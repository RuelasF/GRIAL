# M06: etapas Runge–Kutta con halos CPU

24/oct/2026. Implementación en la copia hidrodinámica `mpi/`. El núcleo ya puede evolucionar bloques radiales mediante MPI en el driver de pruebas. **El programa principal conserva su guarda multiproceso:** sus entradas, observadores, salidas y reinicios distribuidos requieren las etapas posteriores del plan.

## Rutina compartida

`time_integrator.f90` ofrece `advance_rk_step(distributed)`. La llamada sin argumento conserva el camino serial/OpenMP y la selección previa de fronteras CPU/OpenACC. `grial.f90` utiliza esta rutina en lugar de mantener un segundo cuerpo del RK. El argumento lógico es una API interna, no un parámetro nuevo del archivo `.par`.

El llamador calcula dt, atiende la perturbación y actualiza `n_steps` e `integration_time` en el mismo orden previo. La rutina conserva una copia del estado inicial y las tres etapas SSP-RK3:

\[
U^{(1)}=U^n+\Delta t L(P^n),\qquad
U^{(2)}=\tfrac34 U^n+\tfrac14[U^{(1)}+\Delta t L(P^{(1)})],
\]

\[
U^{n+1}=\frac{U^n+2[U^{(2)}+\Delta t L(P^{(2)})]}{3}.
\]

Las expresiones y el orden de operaciones FP64 se trasladaron del programa sin cambios. `cons_to_prim` sigue actuando sobre el estado conservado intermedio, incluidas sus reparaciones. El estado base `u` y el estado intermedio/final `up` mantienen su significado anterior. Tampoco se alteran reconstrucción, Riemann, fuentes, recuperación ni caché métrica.

Con `distributed=.true.` el llamador debe haber inicializado el runtime y `initialize_halo_exchange(layout)` en todos los procesos. Debe existir un reparto común por x/r; los ejes transversales continúan completos. El orden es:

1. Validar CPU, dt común y arreglos locales.
2. Copiar `up` a `u` y llenar primitivas antes del primer RHS.
3. Para cada etapa, calcular RHS, avanzar `up`, recuperar primitivas en celdas físicas y llenar sus halos.
4. Volver al llamador con un RK completo y halos actualizados, antes de diagnósticos, salida o parada.

El llenado utiliza `fill_primitive_halos`: recibe las interfaces internas, aplica fronteras radiales sólo en sus propietarios y después llena aristas y esquinas transversales. Son cuatro llenados por paso completo. Cada llenado usa dos Sendrecv cuando hay varios procesos, con los buffers persistentes de M03. No se llama MPI desde regiones OpenMP.

`parallel_timing:validate_step_size` comprueba dt finito, positivo e idéntico mediante el mismo acuerdo MIN/MAX con marca de error de M05. Los controles de asignación se realizan antes de consultar formas; un estado ausente o de tamaño erróneo en un solo rango rechaza al grupo. Estas validaciones añaden cinco Allreduce por llamada distribuida al RK: CPU, dos para dt, formas y presupuesto radial. No se ha medido ni optimizado escalabilidad.

La ruta serial conserva `accumulate_radial_mass_budget(rk)`, con los pesos previos. La ruta de varios procesos rechaza activar ese presupuesto serial: las caras internas no son fronteras físicas globales y necesitan el tratamiento de D01. Los contadores locales de reconstrucción/Riemann también pueden incluir evaluaciones duplicadas de caras compartidas; no se publican como contadores globales.

## Pruebas de evolución

Desde `~/Codigo/GRIAL/mpi/`:

```bash
make -j2 ARCH_FLAGS='' all test
make -j2 USE_MPI=1 ARCH_FLAGS='' all test test-mpi test-domain-mpi test-halos-mpi test-interfaces-mpi test-time-mpi test-rk-mpi
```

`test-rk` ejecuta el control de un proceso. `test-rk-mpi` ejecuta la matriz local real de 1/2/4 procesos con dos hilos OpenMP por proceso. Los sockets MPI locales del agente requieren salir de su sandbox; no se accede a CADS ni se envían trabajos.

El driver `tests/test_distributed_rk.f90` comienza con SRHD1D y amplía a estos diez fixtures:

| Fixture | Caso | Malla global | Geometría/métrica | Tiempo final |
|---:|---|---|---|---:|
| 1 | Onda de densidad advectada, velocidades transversales no nulas | 73×1×1 | Cartesiana/Minkowski | 0.8 |
| 2 | Sod | 73×1×1 | Cartesiana/Minkowski | 0.15 |
| 3 | Choque fuerte | 73×1×1 | Cartesiana/Minkowski | 0.03 |
| 4 | KHI | 73×1×12 | Cartesiana/Minkowski | 0.03 |
| 5 | Jet | 73×1×12 | Cilíndrica/Minkowski | 0.01 |
| 6 | Michel | 73×1×1 | Esférica logarítmica/EF | 0.2M |
| 7–8 | FD, a=0 y a=0.9 | 73×8×8 | Esferoidal logarítmica/KS | 0.1M |
| 9–10 | FM, a=0 y a=0.9 | 73×16×8 | Esferoidal logarítmica/KS | 0.1M |

Cada fixture combina Godunov, TVD–MC, WENO3, MP5 y WENO5 con HLLE/HLLC: cien trayectorias por tamaño del comunicador. Sensor activo con umbrales por defecto, CFL=0.3, halo de tres capas, pisos rho=1e−10, p=1e−13, D=1e−10 y tau=3e−13. Índice adiabático 5/3 en los tres controles 1D iniciales y 4/3 en los restantes. No se usa MP5 sin sensor como control de choque aceptado.

El reparto de 73 celdas es desigual con 2/4 procesos. La referencia completa se evoluciona con la ruta serial de la misma rutina. El estado físico inicial se copia por índice global desde esa referencia para aislar la evolución del redondeo previo de libm vectorizada al cambiar el largo de las tiras, documentado en M02. No es una nueva prueba de generación de ruido o dato inicial distribuido. Cada fixture de toro debe contener celdas con rho>100 veces el piso antes de evolucionar; un toro vacío hace fallar el test.

Se exigen igualdad exacta de primitivas/conservadas finales y de toda la secuencia dt, tiempo final común, datos finitos, densidad/presión positivas y exactamente tres recuperaciones por celda física y paso. La onda avanza 195 pasos por configuración, con velocidades (0.35,0.2,−0.15); se comprueba conservación global de las cinco variables con tolerancia 1e−11 y ausencia de reparaciones. En total son 2980 pasos completos por proceso y ejecución de cien trayectorias.

Los archivos `rk_cases_rank_*.txt` guardan fixture, esquema/solver, malla local, pasos, los 16 contadores locales, tiempo, errores, mínimos y cinco acumuladores de reparaciones. `rk_rank_*.txt` resume cada rango. El harness comprueba cien recibos por proceso y los conteos de recuperación, sin contar halos como celdas físicas.

La matriz tiene 19 controles: tres ejecuciones normales y ocho errores provocados en 2/4 procesos. Se rechazan dt cero/NaN/distinto, arreglo ausente o de forma errónea, intercambio no inicializado, presupuesto radial serial activado y llamada desde OpenMP. Hay timeout y conservación del log de cada rechazo.

Estos pilotos cortos verifican equivalencia del integrador distribuido. No cierran estabilidad prolongada de toros, crecimiento/saturación PPI, convergencia física de KHI/jet ni validación 3D científica de producción.

## Evidencia y continuidad

`.validation/m06_rk_20261024_v1/` conserva:

- `VERIFIED.json`: suites completas sin MPI/MPI optimizadas y debug con bounds/traps; 29 programas Fortran y dos grupos de integración por build. Matrices MPI runtime 23/dominio 6/halos 23/interfaces 7/tiempo 47/RK 19.
- `HYDRO_VERIFIED.json`: veinte comparaciones del programa completo contra M05 en Sod/KHI/Michel/FD/FM, con referencias de iguales flags y un hilo. Campos, contadores y presupuesto radial idénticos; este control es independiente de comparar dos llamadas a la nueva rutina compartida.
- `FIXTURE_SUPPORT_VERIFIED.json`: los cuatro controles de toro contienen fluido. FD a=0/0.9 tiene 224/320 celdas por encima de 100 veces el piso, y FM a=0/0.9 tiene 208/448. El intento preliminar FM con ocho celdas polares omitía el cuerpo del toro a=0 y se descartó; se conservan fuente/logs/recibos en `coarse_fixture_attempt/` y el rechazo en `fixture_support.log`. Se corrigió sólo el fixture a 16 celdas polares y se añadió la comprobación de soporte. Después de las suites completas se repitieron `test-rk` y las matrices RK en los cuatro builds; los otros módulos y sus controles no cambiaron.
- Controles aceptados, logs, fuentes previas/posteriores, parches y manifiestos. Los incidentes iniciales de compilación del fixture se conservan con sus logs; se corrigieron los nombres de pisos y un intento de asignar la constante del backend, sin cambios al solver.
- Verificación SHA256 de los 52 archivos protegidos del proyecto superior y de la referencia inmutable de 41 archivos.

El programa principal todavía rechaza evolución y `--check` multiproceso antes de leer entradas o escribir archivos. Los drivers no usan sus escritores ni producen un checkpoint distribuido. Las rutas de lectura/configuración/error/E/S del programa deben adaptarse antes de retirar esta guarda; no basta con añadir `mpirun` al comando de una simulación.

El siguiente paso del plan es T01, reparto r×theta manteniendo phi completo. D01–D03 y C01–C04 completarán balances, observadores, perturbación y reinicios/salida antes de campañas MPI largas. R04-CADS conserva su campaña y aceptación independiente; los estados remotos previos mantienen su fecha. No se ejecutó GPU, multinodo, prueba real de 24 h ni benchmark de escalabilidad.

Para revertir este incremento, revisar los cambios posteriores y usar `before/` y `changes.patch` sólo dentro de `mpi/`; el cambio al contexto superior tiene un parche separado. No restaurar el repositorio completo ni modificar `reference/` o campañas congeladas.
