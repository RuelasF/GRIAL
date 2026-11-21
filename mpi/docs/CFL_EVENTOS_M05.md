# M05: CFL global y decisiones temporales

24/oct/2026. Copia hidrodinámica CPU `mpi/`. Se coordinan el CFL, la época y los eventos entre procesos, conservando la política serial. La evolución temporal del fluido con varios procesos sigue bloqueada hasta M06; no hay aún checkpoint ni observadores distribuidos.

## CFL que se conserva

El código utiliza la cota causal de la métrica fija, no una estimación nueva de velocidad del fluido:

\[
\lambda_d=\alpha\sqrt{\gamma^{dd}}+|\beta^d|,\qquad
\Lambda_d=\max_{\mathrm{procesos,celdas}}\lambda_d,
\]

\[
\Delta t=\frac{\mathrm{CFL}}{\Lambda_x/\Delta x+
  \mathbf{1}_{N_y>1}\Lambda_y/\Delta y+
  \mathbf{1}_{N_z>1}\Lambda_z/\Delta z}.
\]

Los espaciados son los lógicos uniformes; el Jacobiano radial logarítmico ya está en la métrica. Se mantiene la contribución x como en la implementación serial, incluso en fixtures degenerados. Para y/z se usan los tamaños globales cuando están disponibles, con el comportamiento local previo para fixtures seriales sin metadatos.

**Se reducen primero los tres máximos direccionales y luego se suma.** Tomar el mínimo de los dt locales puede ser insuficiente si los máximos de distintas direcciones están en bloques distintos. El control sintético no desplazado produce:

| Procesos | dt global 2D | mínimo local 2D | dt global 3D | mínimo local 3D |
|---:|---:|---:|---:|---:|
|1|0.004|0.004|0.0011764705882352942|0.0011764705882352942|
|2|0.004|0.004878048780487805|0.0011764705882352942|0.0015151515151515152|
|4|0.004|0.004878048780487805|0.0011764705882352942|0.0016260162601626016|

Son unidades del control sintético; no constituyen una medida de rendimiento ni una simulación astrofísica.

## Implementación

Nuevo módulo `parallel_timing.f90`, independiente del fluido, de relojes reales y de escritores:

- `validate_cfl_context(spacing,active_axes,courant,cache_ok)`: valida cachés/espaciados y acuerdo del contexto. Se ignoran espaciados inactivos, como antes.
- `reduce_causal_maxima(local_maxima,local_ok,global_maxima)`: reduce con MAX los tres límites y una marca de fallo. Una métrica inválida en un rango impide continuar a todos.
- `cfl_step_from_maxima(...)`: conserva el orden x+y+z y calcula dt, con comprobaciones de tasas/overflow y avance positivo.
- `coordinate_time_actions(time,step,target,event_times,event_enabled,local_requests,decisions)`: verifica época, tiempo final y calendario; acuerda los eventos y solicita checkpoint/parada para el grupo.
- `clip_step_to_target(time,proposal,target)`: conserva la comparación y resta seriales al recortar el último paso; rechaza un paso que no avance el reloj físico.
- `any_time_request(local_request)`: OR colectivo de una petición de checkpoint tomada después de la salida, como en el orden serial.

Las decisiones se devuelven en `time_decisions`: `perturb`, `save`, `checkpoint`, `stop` y `complete`. Los tres eventos son perturbación, salida VTK y checkpoint físico; las dos peticiones locales son parada y checkpoint por tiempo real.

Se exige igualdad exacta de los valores relevantes del calendario, contador de paso y contexto CFL. El contador sigue siendo el entero de la implementación/checkpoint actual y es exactamente representable en FP64 al compararlo. No se comparan calendarios deshabilitados: sus umbrales se normalizan a 0 antes de validar. Un reinicio que ya alcanzó o superó el objetivo termina sin retroceder, como en la ruta previa.

Los metadatos no finitos se sustituyen por valores seguros y una marca de error **antes** de MPI_MIN/MAX. Se comprueba esa marca antes de utilizar el resultado. En la exploración de las celdas se separan los controles de finitud y positividad, porque Fortran no garantiza cortocircuito lógico; así no se compara ni se calcula la raíz de un dato ya identificado como inválido.

El CFL comprueba las formas de lapse, shift e inversa espacial antes del bucle. Los helpers de CFL conservan el uso serial sin inicializar el runtime para las regresiones heredadas; una ejecución MPI debe inicializarlo primero. Las decisiones temporales requieren siempre runtime activo. Todas las llamadas MPI se realizan fuera de regiones OpenMP mediante el guard FUNNELED.

La comparación del contexto utiliza dos reducciones MIN/MAX con marca de fallo y la de velocidades una MAX: tres colectivas por actualización CFL. Las decisiones de cada extremo del RK utilizan tres colectivas; la petición posterior a salida utiliza una. En una iteración completa son diez Allreduce, además de los futuros halos/RK; la salida por objetivo añade su comprobación final. Se prioriza el acuerdo comprobable. No se ha medido ni optimizado escalabilidad en este incremento.

## Conexión al programa

`grial.f90` consulta las decisiones antes de tomar la rama de finalización, de modo que todos participan incluso si alguno tiene una época incoherente. Solicita la parada por el mismo criterio local anterior:

`elapsed_wall + max_observed_full_step >= limit - margin`.

Un OR de solicitudes hace que el grupo pare si cualquiera necesita hacerlo. **No se comparan relojes absolutos de nodos distintos**: cada proceso calcula su petición con sus propios segundos transcurridos y presupuesto. La finalización física tiene prioridad sobre una solicitud de parada en el instante de completar el objetivo.

El orden es:

1. Comprobar época/calendario y petición de parada entre RK completos.
2. Calcular el CFL global y recortar sólo al tiempo final.
3. Aplicar la perturbación si corresponde al inicio del paso, como antes.
4. Completar las tres etapas RK existentes y volver a comprobar época/calendario antes de diagnósticos/salida.
5. Guardar según el evento físico; después medir la petición de checkpoint por tiempo real y propagarla con OR.
6. Conservar el guardado final y el estado `completed`/`stopped_for_walltime` de R01–R04.

No se cambia la cadencia científica: la perturbación se aplica en el primer inicio de paso con t>=t_pert, no se recorta dt para llegar a ella. La salida continúa atendiendo un umbral por paso; no se añaden interpolación temporal ni múltiples salidas por intervalos saltados. El calendario de checkpoints conserva la actualización serial, también cuando el guardado lo solicita el reloj real.

La medida del paso mantiene el instante tomado al comienzo, por lo que incluye el tiempo de coordinación y salida al reservar el paso completo más lento. Estas decisiones no sustituyen la calibración real del margen de R04-CADS/V04.

No hay nuevos flags de usuario, señales, requeue ni envíos Slurm. Los cuerpos de recuperación, flujos, reconstrucción y RK no cambian. El programa todavía rechaza varios procesos antes de leer/evolucionar el problema: sólo los drivers de prueba utilizan estas capacidades MPI. Los escritores de un proceso no se habilitan para MPI mediante este incremento.

## Pruebas

Desde `~/Codigo/GRIAL/mpi/`:

```bash
make -j2 ARCH_FLAGS='' all test
make -j2 USE_MPI=1 ARCH_FLAGS='' all test test-mpi test-domain-mpi test-halos-mpi test-interfaces-mpi test-time-mpi
```

`test-time` ejecuta los controles de un proceso y `test-time-mpi` lanza 1/2/4 procesos reales locales con dos hilos OpenMP. El sandbox de un agente puede impedir los sockets locales; las pruebas se ejecutaron fuera de él cuando fue necesario, sin acceder al clúster.

`test_parallel_timing.f90` verifica:

- Seis controles sintéticos 1D/2D/3D, con y sin lapse/shift, colocando máximos direccionales en bloques diferentes. Compara con la fórmula serial y demuestra el fallo del mínimo de dt locales.
- Cinco geometrías con caché completa y por bloques: Minkowski cartesiana, cilíndrica 2D, EF radial, KS logarítmica meridional y KS logarítmica 3D con spin 0.9. El dt coincide exactamente al repartir N=37 en 1/2/4 procesos.
- Dos calendarios simulados de cuatro pasos, sin fluido: propuestas 0.25 y 0.3 hasta 1, perturbación en 0.375 atendida una sola vez en 0.5/0.6, salidas y checkpoints con la cadencia anterior, y un checkpoint solicitado por un solo rango.
- Petición de parada/checkpoint de un único rango, prioridad de finalización física, eventos deshabilitados, reinicio ya completo y ausencia de solicitudes espurias.
- Fallos en un rango en 2/4 procesos: espaciados, CFL o direcciones distintos; NaN, lapse/inversa inválidos o shift infinito; cachés ausentes o formas erróneas; época/paso/objetivo/calendario/flags distintos; datos temporales no finitos, paso sin avance, tasa no representable y llamada desde OpenMP. La matriz tiene 47 controles, timeout y recibos completos, sin datos de simulación ni checkpoint MPI.

Los NaN y rechazos se probaron también con traps invalid/zero/overflow. No se afirma validación de cualquier métrica o rango numérico posible a partir de esta matriz acotada.

## Evidencia y límites

`.validation/m05_time_20261024_v1/` conserva:

- `VERIFIED.json`: 28 regresiones Fortran y dos grupos de integración en cuatro builds sin MPI/MPI optimizado/debug. Matrices MPI 23/6/23/7/47 para infraestructura/dominio/halos/interfaces/tiempo.
- `CFL_DEFAULT_VERIFIED.json`: once dt seriales de controles sintéticos y cachés reales, idénticos porbytes frente a `before/evolution.f90`; también comprueba el uso heredado sin `parallel_initialize`.
- `HYDRO_VERIFIED.json`: 20 comparaciones cortas de evolución de un proceso en Sod/KHI/Michel/FD/FM con referencia pre-M05 de iguales flags, exactas en campos, contadores y presupuestos. Permanecen separadas las diferencias O0/O3 y sumas atómicas OpenMP ya documentadas.
- Fuentes previas, parche, manifiestos de fuentes/ejecutables y logs. Referencia inmutable y archivos protegidos del proyecto superior verificados porSHA256.

Los calendarios y dt comunes no implementan comunicación de las etapas RK, perturbaciones aleatorias globales, contadores/observadores globales ni checkpoint distribuido. Tampoco validan una ejecución física multinodo, GPU o continuidad real de 24 h. La campaña R04-CADS conserva fuentes, trabajos y aceptación independiente; los estados remotos anteriores mantienen fecha.

**Siguiente M06:** conectar halos y reloj común en cada etapa, comprobar equivalencia de evoluciones SRHD1D antes de ampliar a 2D/GRHD y revisar las rutas fatales/configuración y salida antes de retirar la guarda del ejecutable. No se cierra Control M por haber cerrado sólo M05.

Para revertir, revisar cambios posteriores y usar `before/` y `changes.patch` sólo sobre este incremento en `mpi/`. No restaurar globalmente ni editar `reference/` o fuentes congeladas de CADS.
