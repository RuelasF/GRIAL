# D01 — diagnósticos globales y presupuesto conservativo

**Actualizaciones posteriores D02/D03:** el cierre y las cifras D01 siguientes son históricos. La versión actual admite `root_only=.true.` en los dos colectores y lo utiliza para publicar; conserva Allreduce de la escala L2 y los resultados globales por defecto. PPI/acreción/ondas cuentan con [su cierre local D02](OBSERVABLES_D02.md). [D03](PERTURBACIONES_D03.md) añade el cambio instantáneo de perturbación y amplía el presupuesto a 62 columnas; reinicios y producción MPI siguen pendientes.

Proyecto interno hidrodinámico CPU de GRIAL, cierre original del 24/oct/2026. Implementación en `parallel_diagnostics.f90`, con captura en `evolution.f90` y `time_integrator.f90`, y publicación de series desde rango 0 en `output.f90`. Evidencia original en `.validation/d01_diagnostics_20261024_v1/VERIFIED.json`. El programa principal conserva su guarda multiproceso; las extensiones posteriores se identifican arriba.

## Cantidades y unidades

Se integran las cinco variables evolucionadas, en el orden **D, tau, S1, S2, S3**. Son densidades respecto a la malla lógica: el Jacobiano y sqrt(gamma) ya están incluidos. Cada celda física contribuye una vez, sin halos, multiplicando por `dx*dy*dz`. No se aplica otra vez sqrt(gamma).

`totals(eq_de)` es la masa bariónica del dominio computacional. En las coordenadas esferoidales de los toros, `totals(eq_vz)` es J, la integral de S_phi. En cartesianas S3 es el momento covariante de la tercera dirección, no un momento angular. Los totales incluyen el dominio completo, también su parte interior al horizonte si está en la malla. No equivalen a una tasa de acreción en la sonda del horizonte, calculada por separado en D02.

El estado global contiene además masa por encima de `1000*rho_floor`, extremos de rho/p, máximo de gamma_ij v^i v^j, celdas de atmósfera/invalidación, normas de U y contadores int64. La clasificación de estado inválido considera primitivas no positivas y velocidad fuera del rango físico; los datos no finitos se rechazan colectivamente antes de reducirlos.

Las normas son medias sobre las celdas de la malla lógica: L1=mean(abs(U)), L2=sqrt(mean(U**2)), Linf=max(abs(U)). El cálculo L2 se escala por el máximo global antes de elevar al cuadrado. `conserved_error_norms(reference,...)` aplica las mismas definiciones a U-U_ref, con una referencia local de igual forma. No son órdenes de convergencia ni normas respecto a una solución exacta si no se proporciona esa referencia.

## Balance del paso RK

Para Q=sum(U*DeltaV), la comprobación es

    R = Q - Q_initial - sum(I_F_lower - I_F_upper) - I_source - I_recovery - I_perturbation.

Desde D03, `I_perturbation` registra por separado `sum((U_post-U_pre)*DeltaV)` al inyectar, fuera del RK y sin pesos temporales. Es cero en los controles D01 originales. El cambio impuesto al fluido no se atribuye al recuperador ni a las fuentes geométricas.

Los flujos son los flujos numéricos que realmente utiliza el RHS, con área lógica transversal. Sólo se incorporan caras externas no periódicas de cada eje y únicamente desde sus propietarios. Las interfaces MPI se excluyen. Los extremos periódicos no se contabilizan como flujo externo. En los polos se conserva el flujo de área nula del solver existente.

Las fuentes se capturan de la misma evaluación `calc_sources` utilizada por el RHS. Su suma local usa una reducción OpenMP, sin comunicaciones dentro de la región paralela. Las tasas de flujo y fuente se acumulan con dt y los pesos SSP-RK3:

    RHS: [1/6, 1/6, 2/3].

Una reparación de conservadas después de recuperar primitivas modifica la etapa ya actualizada. Su contribución al estado final tiene otros pesos:

    recovery: [1/6, 2/3, 1].

Se obtiene el cambio de cada celda inmediatamente antes/después de `cons_to_prim`, y luego se reduce localmente con OpenMP. Esto evita inferir diferencias pequeñas restando acumuladores históricos grandes. El diagnóstico `raw_repair` conserva la suma histórica sin ponderar; no debe sustituir a `I_recovery` en el balance del paso completo.

Se rechaza consultar/publicar un presupuesto a mitad de RK. Al cerrar el paso se exige que se hayan capturado las tres evaluaciones de RHS y sus tres recuperaciones, en orden. Este módulo registra el paso; el llamador sigue manteniendo `integration_time` y `n_steps`.

## Contadores de trabajo y de celdas

Se suman los contadores con un tipo MPI entero de ocho bytes obtenido mediante `MPI_Type_match_size`. El caso sintético utiliza valores mayores que 2^31 por rango y comprueba la suma exacta.

Los eventos C2P 8:13 pertenecen a recuperaciones de celdas físicas. Las reconstrucciones, los estados de cara y las llamadas HLLC son evaluaciones de trabajo; incluyen caras compartidas calculadas en ambos bloques. No deben interpretarse como celdas únicas ni como reparaciones adicionales del fluido.

Se publican por separado:

- `unique_faces_per_rhs`: posiciones de cara globales, incluyendo ambos índices extremos como en el cálculo serial.
- `work_faces_per_rhs`: posiciones evaluadas por todos los bloques.
- `extra_interface_faces`: evaluaciones adicionales por interfaces MPI desde la inicialización del presupuesto, incluyendo tres etapas por paso.

En cada eje activo, el primer recuento es `(N_d+1)*product(N_transverse)`. Al repartir un eje entre P_d bloques, el trabajo usa `(N_d+P_d)*product(N_transverse)`. El tercer eje permanece completo. Este recuento no intenta convertir las llamadas HLLC en caras únicas: en los polos y ejes degenerados algunas caras no invocan ese solver.

## API y uso

| Rutina o tipo | Función |
| --- | --- |
| `state_diagnostics` | Estado global, normas, extremos, contadores y trabajo de caras |
| `conservation_balance` | Q inicial/actual, flujos por variable/lado/eje, fuentes, reparaciones, perturbación y residuo |
| `initialize_parallel_diagnostics(layout)` | Valida configuración y guarda la línea base del presupuesto; activa los hooks CPU |
| `collect_global_state(state)` | Obtiene el estado global en todos los rangos |
| `conserved_error_norms(reference,l1,l2,linf)` | Obtiene errores globales respecto a una referencia local |
| `collect_conservation_balance(balance)` | Reduce acumuladores y forma el residuo global |
| `accumulate_perturbation_budget(change)` | Acumula el cambio local de U por perturbación, fuera del RK; multiplica por el volumen lógico |
| `finalize_parallel_diagnostics()` | Desactiva el presupuesto al acabar un paso completo |
| `write_distributed_diagnostics(time)` | Llamada colectiva: recoge datos y escribe únicamente desde rango 0 |

Los hooks `begin_budget_step`, `reset_rhs_budget`, `capture_boundary_fluxes`, `capture_source_sum`, `accumulate_rhs_budget`, `accumulate_recovery_budget` y `end_budget_step` los coordinan el integrador/RHS compartidos. No se debe llamar al RHS de manera aislada con el presupuesto activo, fuera de esa secuencia. La captura de caras puede ejecutarse dentro del bucle OpenMP; no realiza llamadas MPI.

Inicializar después de construir la malla, caché y estado físico. Hacen falta `do_numerical_diagnostics=true`, CPU y runtime inicializado. El inicializador no borra los contadores del llamador: en los controles se ponen a cero después de preparar el dato inicial, excluyendo el trabajo de arranque. Los acumuladores del presupuesto y el contador de trabajo adicional empiezan con esta inicialización.

La topología deriva de las fronteras de cada caso. Conv/KHI son periódicos; los choques y Michel tienen extremos radiales externos; FD/FM sagitales/3D tienen fronteras polares y phi periódico; el jet tiene fronteras externas, sin imponer periodicidad Y. El soporte del intercambio polar sigue limitado por T01, independientemente de que el colector pueda comprobar algebraicamente una partición.

## Archivos de salida

Bajo `output_folder/scheme_name/`:

| Archivo | Contenido |
| --- | --- |
| `mpi_global_diagnostics.dat` | Masa, masa del cuerpo, S3, extremos, normas de U y recuentos de celdas |
| `mpi_conservation_budget.dat` | Q inicial/actual, integral de flujo `(variable,lado,eje)`, fuentes, reparación ponderada, perturbación y residuo; 62 columnas desde D03 |
| `mpi_numerical_diagnostics.dat` | Contadores globales, reparación cruda y recuentos de trabajo en caras |

La carpeta debe existir y todos los rangos deben llamar a la publicación. Sólo rango 0 abre los archivos. La primera muestra exige archivos nuevos; las siguientes exigen tiempo y paso crecientes. Una escritura fallida o duplicada produce un error colectivo y conserva la evidencia disponible. No hay protocolo transaccional de reinicio para estas tres series: `do_restart=true` se rechaza al inicializar D01 hasta C01/C03. Checkpoints V9 y series seriales anteriores conservan su formato y comportamiento por defecto.

## Controles y alcance

Ocho repartos: 1×1,2×1,1×2,4×1,1×4,2×2,4×2,2×4, dos hilos por rango. Los problemas con ny=1 reparten sólo radio entre todos los procesos. Configuraciones MP5/WENO5 con sensor y HLLE/HLLC:

| Caso | Malla global | Tiempo final |
| --- | --- | --- |
| Onda periódica 3D | 37×19×8 | 0.01 |
| Sod, Gamma=1.4 | 73×1×1 | 0.4 |
| Choque fuerte, Gamma=5/3 | 73×1×1 | 0.4 |
| Michel | 73×1×1 | 0.2M |
| FD 3D, a=0 | 37×19×8 | 0.03M |
| FM 3D, a=0.9 | 37×19×8 | 0.03M |
| Reparación forzada | 37×19×8 | 0.0001 |
| Jet axisimétrico cilíndrico | 37×1×19 | 0.01 |

Se añaden dos choques fuertes Godunov N400 hasta 0.4, uno con cada solver: **34 trayectorias por grid/build, 1936 pasos por rango sumados entre esos casos, 544 comparaciones opt/debug**. Los estados iniciales se copian por índices de un control completo para aislar el presupuesto/MPI. El caso de reparación forzada utiliza deliberadamente conservadas incompatibles con primitivas de atmósfera: es un control del peso de las correcciones, no un dato inicial físico de producción.

También se prueban sumas/extremos/normas sintéticos, halos envenenados excluidos del diagnóstico, int64, máscaras externas en los tres ejes y signos de flujo en Conv/Sod/Jet. Quince rechazos en 2/4 procesos: metadatos/malla, geometría, pisos, contadores desactivados, reinicio, estado sin inicializar, etapas/época/dt, presupuesto parcial, referencia mal dimensionada, escritura duplicada, no finitos y rango de contadores. El postanálisis vuelve a auditar los mismos archivos para exigir por separado masa y J, sin repetir evoluciones.

Criterio 1e-11 conservado. En optimizado/debug, máximo residuo relativo vectorial 3.002e-14, masa 2.817e-14 y J en Kerr 2.238e-15; campos y secuencia dt exactamente iguales al control completo. Los Godunov N400 cambian la masa relativa alrededor de 1.1e-11; el presupuesto con flujo numérico real deja residuos de masa alrededor de 1e-16. Se distingue el flujo de frontera de una pérdida numérica o una reparación.

Las suites completas incluyen 31 programas Fortran y dos grupos de integración por build sin MPI/MPI, optimizado/debug. Las matrices MPI anteriores y 28 controles del programa habitual se conservan en la evidencia. Los intentos de enlace de pruebas, el identificador incorrecto `strong` en un archivo de control y la revisión de máscaras del jet se archivan con sus correcciones.

## Reproducción

```bash
make -j2 ARCH_FLAGS='' test-diagnostics
make -j1 USE_MPI=1 ARCH_FLAGS='' test-diagnostics-mpi
```

Ejecutar matrices MPI una por vez; sus controles incluyen el coste de recalcular la referencia completa en cada rango y no son benchmarks. El harness crea carpetas distintas para cada lanzamiento, conserva salidas y verifica las series escritas desde rango 0. El dossier de evidencia incluye fuentes antes/después, parches, ejecutables y manifiestos SHA256.

D01 cierra el colector/presupuesto CPU local. La extracción PPI, acreción y ondas distribuida pertenece a D02; perturbaciones a D03; checkpoints/redistribución/continuidad y visualización a C01–C04. La entrada multiproceso de producción sigue bloqueada. No se ha ejecutado GPU, multinodo ni una prueba de escalabilidad o continuidad real de 24 horas en este incremento.
