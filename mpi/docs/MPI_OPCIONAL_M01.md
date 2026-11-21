# M01 — Compilación MPI opcional e infraestructura de ejecución

Implementado localmente el 24/oct/2026 UTC, por autorización de avanzar mientras la campaña de reinicio CADS conserva su propia versión congelada. El código principal con SRMHD y los archivos de CADS permanecen fuera de este incremento.

## Compilar y ejecutar

Desde `~/Codigo/GRIAL/mpi/`:

```bash
# Ruta habitual: un proceso, OpenMP, sin enlazar MPI.
make -j2 ARCH_FLAGS=''
OMP_NUM_THREADS=2 ./grial archivo.par

# Ruta opcional: requiere un wrapper Fortran con mpi_f08.
make -j2 USE_MPI=1 ARCH_FLAGS=''
OMP_NUM_THREADS=2 mpirun -n 1 ./grial_mpi archivo.par

# Arranque/cierre y topología, sin parámetros ni salidas científicas.
OMP_NUM_THREADS=2 mpirun -n 4 ./grial_mpi --parallel-check

# Regresiones y controles de la infraestructura.
make test ARCH_FLAGS=''
make USE_MPI=1 test test-mpi ARCH_FLAGS=''
```

`MPIFC` selecciona el wrapper MPI (por defecto `mpifort`). `PYTHON` selecciona el intérprete de los controles. Compilación/ejecución verificadas con GNU Fortran y Open MPI 5.0.7 locales. La disponibilidad/compatibilidad del entorno MPI de CADS/TAMU no se infiere de esta prueba.

El ejecutable normal sigue siendo `grial`. La variante MPI está en `.build/mpi/grial_mpi`, expuesta mediante `grial_mpi`; sus objetos, módulos y binarios de prueba se construyen desde otro directorio. Esto evita leer accidentalmente los `.mod` del build sin MPI. Los modos pueden alternarse sin limpiar el otro. Hay que usar copias de build separadas o limpiar el modo correspondiente si se cambian compilador/flags, como antes.

`USE_MPI` acepta 0/1. MPI+OpenACC se rechaza en esta etapa CPU. La variante sin MPI no requiere `mpifort`, `mpi_f08` ni bibliotecas MPI: el enlace se comprobó con `ldd`.

## Alcance efectivo

**La evolución todavía usa una malla completa en un solo proceso.** M01 no implementa descomposición, halos, CFL global, diagnósticos distribuidos ni E/S MPI. No representa una aceleración distribuida.

- `mpirun -n 1 ./grial_mpi archivo.par` y `--check archivo.par` son válidos.
- Con varios procesos se admite `--parallel-check` y el programa de regresión del runtime.
- Intentar evolucionar o validar un caso con varios procesos aborta antes de leer parámetros, alojar la malla o crear archivos científicos. El MPI binario no ejecuta copias completas independientes de la simulación.
- El binario sin MPI también rechaza lanzamientos que indiquen varios procesos mediante `OMPI_COMM_WORLD_SIZE` o `PMI_SIZE`. No debe usarse un binario sin MPI como programa multiproceso bajo otros lanzadores.

La partición de la malla corresponde a M02; los halos a M03; la evolución distribuida se habilitará al verificar M06. El guard de esta etapa se mantiene hasta entonces.

## API y errores

Nuevo módulo `parallel_runtime.f90`, aislado de los módulos físicos:

|API|Contrato|
|---|---|
|`parallel_initialize()`|Solicita `MPI_THREAD_FUNNELED`, comprueba el nivel recibido y obtiene rango/tamaño; inicialización única propiedad de GRIAL.|
|`parallel_finalize()`|Cierra normalmente el runtime después de la limpieza y en las rutas `--check`/`--parallel-check`.|
|`parallel_require_collective(ok,message)`|Todos los rangos participan en validación mediante `MPI_Allreduce`/`MPI_LAND`; si uno falla se aborta el grupo.|
|`parallel_check_mpi(status,operation)`|Comprueba el estado devuelto por MPI, obtiene su descripción y aborta ante error.|
|`parallel_abort(message,code)`|Fallo desde un rango; `MPI_Abort` termina el grupo, incluidos los compañeros que esperan. Sin MPI termina con error.|

Rango, tamaño y nivel de hilos se exportan como estado protegido, con valores 0/1 en el fallback sin MPI. Las comunicaciones usan `mpi_f08` y `MPI_COMM_WORLD`, con `MPI_ERRORS_RETURN`. No se intenta recuperar una ejecución después de un error MPI.

Las llamadas se realizan desde el hilo principal, fuera de regiones OpenMP activas. Se rechaza una llamada desde una región activa antes de entrar a MPI; el error local resultante también provoca salida del job en el lanzador probado. La validación colectiva requiere que todos sus participantes lleguen al mismo punto; no sustituye un aborto inmediato cuando hay un fallo local que impide llegar.

Los errores `stop/error stop` históricos dentro de módulos físicos/lectores no se reescribieron masivamente. Sólo son alcanzables con un proceso en esta etapa. Deben revisarse antes de habilitar la evolución distribuida para impedir esperas asimétricas. Los fallos anteriores a obtener un contexto MPI válido dependen del manejo inicial de la biblioteca/lanzador.

## Entorno de punto flotante

La instalación local de Open MPI puede levantar excepciones dentro de `MPI_Init_thread`; con traps activos se reprodujo SIGFPE antes de entrar al código de evolución. Se conserva el estado IEEE y se desactiva halting sólo alrededor de inicialización/finalización de la biblioteca; inmediatamente después se restaura el estado original.

Los traps del código hidrodinámico siguen activos. Las pruebas verifican igualdad de flags y halting antes/después y su conservación en cada hilo OpenMP. No se modifica FP64, redondeo ni algoritmos de recuperación/flujo para acomodar MPI.

Open MPI 5.0.7 tiene un comportamiento específico: si `OMPI_MPI_THREAD_LEVEL` difiere del nivel solicitado, su inicialización puede abortar con `MPI_ERR_ARG`, incluso para un nivel superior. Los controles usan la negociación normal; comprueban el rechazo del override0 y la aceptación del override1. No se debe fijar un override distinto para esta versión.

Referencias primarias: [MPI_Init_thread y FUNNELED](https://docs.open-mpi.org/en/v5.0.0/man-openmpi/man3/MPI_Init_thread.3.html), [implementación de Init_thread en5.0.7](https://github.com/open-mpi/ompi/blob/v5.0.7/ompi/mpi/c/init_thread.c), [estado IEEE de Fortran](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2026-1/ieee-set-status.html).

## Verificación y límites de las comparaciones

Evidencia: `.validation/m01_optional_mpi_20261024_v1/VERIFIED.json`, fuentes previas, referencia ejecutable R04 y su build debug, parches y logs.

- 24 programas Fortran y dos grupos de integración en cada modo sinMPI/MPI y optimizado/debug. Se conservan las23 regresiones hidrodinámicas anteriores y se añade el runtime.
- Matriz MPI optimizada/debug con 1/2/4 procesos: vida del runtime, dos hilos OpenMP por rango, error colectivo, aborto local con compañeros esperando, error real de MPI por destino inválido, llamada desde worker, CLI incorrecta, niveles de hilo y rechazo de evolución/`--check` con varios procesos. También rechazo de lanzar el binario sin MPI con 2/4 procesos. Los controles tienen timeout y limpieza de su propio grupo de procesos para detectar bloqueos.
- Sod, Michel, KHI y FD/FM3D comparados contra el código anterior:25 comparaciones con un hilo, usando referencias O3/O0 equivalentes, con campos, tiempo/paso, contadores y presupuestos exactamente iguales. Otras 25 con dos hilos conservan exactamente los campos y contadores enteros; las sumas flotantes de reparación se registran por separado.

Se identificaron dos efectos anteriores a M01, sin relajar tolerancias:

1. KHI tiene diferencia O0/O3 de 3.553e−15 ya en la referencia previa. Por eso se compara cada build con su referencia de las mismas flags.
2. Las sumas atómicas de reparaciones con varios hilos dependen del orden de suma. Dos ejecuciones del mismo binario previo conservaron campos y contadores enteros, con diferencias de hasta 3.553e−14 en esas sumas. La aceptación exacta de acumuladores se comprueba con un hilo; el control de campos se repite con dos. No se cambia el algoritmo de suma en esta tarea: su tratamiento distribuido corresponde a D01.

Estos controles verifican la incorporación del runtime. No sustituyen convergencia física, validación3D prolongada o escalabilidad. No se han usado créditos TAMU ni enviado nuevos jobs para M01. La campaña R04-CADS sigue teniendo su versión y estado propios; su aceptación real permanece pendiente hasta revisar sus recibos/análisis.

**Siguiente: M02**, tamaños/índices/coordenadas globales y bloques radiales locales, manteniendo la evolución multiproceso bloqueada hasta verificar el resto de la ruta.
