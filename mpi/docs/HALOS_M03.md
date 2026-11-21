# M03: halos de primitivas y fronteras físicas

24/oct/2026. Implementación CPU en la copia hidrodinámica `mpi/`. Se intercambian primitivas FP64 entre bloques del primer eje x/r. La evolución distribuida sigue bloqueada en el ejecutable hasta M06.

## Qué se implementó

`halo_exchange.f90` ofrece cuatro rutinas:

- `initialize_halo_exchange(layout)`: comprueba colectivamente el dominio y prepara los buffers.
- `exchange_primitive_halos(p)`: recibe las capas x/r, sin aplicar fronteras físicas ni completar las direcciones transversales.
- `fill_primitive_halos(p)`: intercambia y después llama a las fronteras físicas con los propietarios correctos.
- `finalize_halo_exchange()`: libera los buffers; debe llamarse antes de finalizar el runtime paralelo.

Se conserva `[x,y,z,variable]`, `real*8` y las cinco componentes hidrodinámicas. No hay comunicación de campos magnéticos ni GPU distribuida.

## Preparación y comunicación

La preparación requiere el runtime activo, la malla local de M02, la caché métrica cuando la frontera la utiliza y el inyector de Michel únicamente en el proceso que posee el extremo exterior. Comprueba:

- Layout canónico, rango/tamaño de `MPI_COMM_WORLD`, profundidad de halo y tamaños/offset globales y locales.
- Mismo tamaño global, profundidad, número de variables, periodicidad y nombre de frontera entre procesos.
- Periodicidad radial para KHI/Conv; extremos físicos para los demás casos soportados.
- Caché local, inyector exterior y phi par para los casos con cruce polar.
- Conteo MPI representable sin desbordamiento y resultado de la reserva de memoria.

Sólo después se reservan dos buffers persistentes de `nghost*ny*nz*neq` elementos cada uno. La prueba del conteo inválido rechaza el problema antes de intentar reservar un arreglo enorme. La interfaz usa conteos MPI enteros convencionales; una capa que los exceda se rechaza, no se fragmenta automáticamente.

Cada intercambio empaqueta únicamente celdas físicas transversales `1:ny,1:nz`. Hace dos `MPI_Sendrecv` bloqueantes sobre `MPI_COMM_WORLD`, con tags 201/202:

1. Envía las primeras capas al vecino inferior y recibe las primeras capas del vecino superior en el halo superior.
2. Envía las últimas capas al vecino superior y recibe las últimas capas del vecino inferior en el halo inferior.

La selección mantiene el orden creciente del índice global. Los vecinos inexistentes del dominio (-1) se traducen a `MPI_PROC_NULL`; no se escribe el halo correspondiente al recibir de ese extremo. Se revisa el código de cada llamada y el número de elementos recibidos. Los buffers se reutilizan; no se reservan de nuevo en cada etapa.

No hay reducciones globales en un intercambio normal. Se comprueban localmente el contexto, la forma del arreglo y la continuidad de la configuración; un fallo llama a `MPI_Abort` y termina el grupo. `parallel_require_active` también impide llamar desde una región OpenMP antes de llegar a MPI. Las reducciones de preparación no sustituyen el futuro CFL ni el control temporal global de M05/M06.

Con un único proceso, la periodicidad se resuelve por copia local usando módulo del índice físico. Esto admite `nx=1` y profundidad3 sin leer capas todavía no inicializadas. Las fronteras no periódicas quedan para su propietario. Con varios procesos se mantiene el mínimo de bloque de M02.

## Orden de fronteras y esquinas

El orden del wrapper es **intercambio radial → fronteras físicas radiales de sus propietarios → fronteras transversales**. Los datos recién recibidos se propagan a aristas y esquinas al completar las direcciones transversales. En el reparto radial actual, theta y phi están completos en cada proceso y se conservan las reglas de reflexión/cruce polar y periodicidad existentes.

`set_boundary_conditions` incorpora dos argumentos opcionales:

```fortran
call set_boundary_conditions(p,fill_lower_x,fill_upper_x)
```

Por omisión ambos valen verdadero y se conserva el comportamiento serial. El wrapper pasa `layout%physical_lower(1)` y `%physical_upper(1)`. Esto evita aplicar outflow, reflexión, inyector o periodicidad local sobre una interfaz MPI.

El jet requiere además que su reinicialización de fantasmas a ambiente distinga un extremo físico de una interfaz MPI: conserva las primitivas recibidas en x con índices transversales físicos. Luego aplica las fronteras transversales con su precedencia original. En Michel, sólo el propietario exterior necesita acceder al inyector; los demás procesos no lo reservan en la prueba.

Las capas inferiores útiles son `1-nghost:0`; con profundidad3 son `-2,-1,0`. Las superiores son `nx+1:nx+nghost`. El índice adicional `-nghost` reservado por la implementación histórica no se transmite. Las fronteras transversales pueden escribir sobre ese índice adicional como antes; no debe confundirse con una cuarta capa comunicada.

**No llamar a las fronteras seriales sin máscaras después del intercambio MPI:** sobrescribiría las interfaces. La conexión al integrador sustituirá sus llamadas por el wrapper en cada etapa, una vez verificados M04 y M05. Actualmente el programa sigue usando la ruta de un proceso.

## Pruebas y comandos

Desde `~/Codigo/GRIAL/mpi/`:

```bash
make -j2 ARCH_FLAGS='' all test
make -j2 USE_MPI=1 ARCH_FLAGS='' all test test-mpi test-domain-mpi test-halos-mpi
```

`test-halos` comprueba el caso de un proceso. `test-halos-mpi` ejecuta procesos reales locales 1/2/4 con dos hilos OpenMP por proceso; no envía trabajos Slurm. Los sockets locales requieren ejecución fuera del sandbox de un agente cuando éste los restringe.

Se prueban los 15 tipos de frontera, tanto por nombre como por identificador numérico, con profundidades 1/2/3 y N=37 (reparto desigual). Son 90 configuraciones por proceso, más tres controles degenerados `nx=1` con un proceso. Cada configuración usa dos conjuntos sucesivos de datos codificados por índice global y componente:

- La prueba del intercambio puro verifica cada entrada asignada, incluidos los valores que deben permanecer intactos. Detecta inversión de capas, componentes cruzadas y comunicación accidental de fantasmas transversales o del índice extra.
- El wrapper se compara exactamente con el corte correspondiente de la malla serial después de aplicar sus fronteras, incluidos extremos físicos, aristas y esquinas.
- Se incluyen presión/densidad bajo pisos en posiciones codificadas, velocidades de ambos signos, geometría cartesiana/cilíndrica, EF/KS, Michel, polvo, FD/FM y cruces polares. Son controles de transporte/fronteras, no evoluciones físicas de esos estados sintéticos.
- Se vuelve a llenar la malla con otros datos y se repite para detectar buffers viejos.

La matriz incluye rechazos en 2/4 procesos por layout inválido, periodicidad incompatible, formas globales distintas, arreglo local incorrecto, llamada desde OpenMP, contexto no inicializado, conteo MPI excesivo, caché o inyector ausentes y nz impar en el cruce polar. Cada prueba tiene timeout y termina su propio grupo si se bloquea; no se considera válida una ejecución que sólo termine sin producir su recibo completo.

La captura por proceso usa archivos pequeños de prueba: el encaminamiento stdout de MPI puede intercalar líneas largas. No se crean campos de simulación ni checkpoints en esa matriz.

## Evidencia y alcance

`.validation/m03_halos_20261024_v1/` conserva fuentes previas, parche, manifiestos y logs. `VERIFIED.json` registra el cierre comprobado.

- 26 regresiones Fortran+dos grupos de integración en cada build sinMPI/MPI optimizado/debug. Debug con límites y traps de invalid/zero/overflow.
- Matrices M01 de 23 controles, M02 de 6 y M03 de 23 en cada buildMPI; halos y fronteras comparados con error0.
- 192 estados de fronteras por defecto —también los índices extra asignados— comparados porbytes con una compilación de las fuentes pre-M03, sin el módulo de halos. Se conservan los dos volcados y suSHA256; no se utiliza la propia frontera modificada como única referencia.
- 20 comparaciones cortas de evolución de un proceso (Sod/KHI/Michel/FD/FM) con referencias de iguales flags y un hilo, incluyendo contadores y presupuesto radial. Las diferencias atómicas OpenMP y O0/O3 ya documentadas permanecen fuera de esta exigencia de identidad porbytes.

No se modifican reconstrucción, Riemann, fuentes, recuperación de primitivas ni RK para acomodar MPI. La referencia inmutable y el código principal se verifican por SHA256. Las campañas CADS mantienen sus fuentes, trabajos y aceptación independiente; sus estados anteriores conservan fecha.

## Siguiente paso y reversión

**M04:** hacer que las protecciones de reconstrucción y el sensor distingan bordes físicos de interfaces MPI y verificar estados/flujos compartidos. Después M05 CFL/decisiones globales y M06 RK distribuido. Primero se cerrará SRHD1D, después SRHD2D y finalmente evoluciónGRHD. Esta etapa no valida una evolución multinodo ni mide escalabilidad.

Para revertir, revisar modificaciones posteriores y usar `before/` y `changes.patch` sólo sobre este incremento en `mpi/`. No restaurar globalmente el proyecto ni modificar `reference/` o fuentes congeladas de CADS.
