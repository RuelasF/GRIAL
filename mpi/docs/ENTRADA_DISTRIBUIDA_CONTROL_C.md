# Control C: entrada distribuida desde GRIAL

Implementación integrada y verificada localmente el 27/oct/2026. La aceptación final se registra en [VERIFIED.json](../.validation/control_c_entry_20261027_v1/VERIFIED.json); no inferirla de la existencia del ejecutable. Esta ruta pertenece sólo a la copia hidrodinámica `mpi/`.

## Compilar y ejecutar

```bash
cd ~/Codigo/GRIAL/mpi
make -j2 USE_MPI=1 ARCH_FLAGS=''
OMP_NUM_THREADS=2 mpirun -n 4 ./grial_mpi --distributed caso.par --grid 2 2
```

Se utiliza el `.par` habitual. En esta etapa la ruta requiere CPU/hidrodinámica, `numerical_diagnostics=true`, strides de diagnóstico/extracción positivos y presupuesto radial serial desactivado. Esas condiciones se comprueban antes de escribir. Los controles físicos y métodos son los existentes; no se selecciona otro solver al activar esta entrada.

Los rangos deben compartir el mismo directorio de trabajo y árbol de archivos, también al recuperar historias y piezas. El acuerdo de argumentos/contenido no comprueba por sí mismo que rutas relativas desde carpetas distintas apunten al mismo árbol compartido.

`--grid nr ntheta` determina el reparto radial/polar; su producto debe coincidir con los procesos MPI y los bloques deben admitir el stencil. Sin esa opción se usa `[nprocesos,1]`. Phi permanece completo en cada bloque. En SRHD 1D y KHI/Jet x–z se usa, por ejemplo, `--grid 4 1`; `--grid 2 2` corresponde a una malla con theta activa. `OMP_NUM_THREADS` son hilos por proceso: cuatro procesos con dos hilos requieren ocho cores para evitar sobreasignación. Las pruebas locales con oversubscription no representan rendimiento de producción.

La compilación sin MPI también permite `./grial --distributed caso.par` con un solo proceso para verificar el mismo camino de cálculo/salida. Un binario sin MPI continúa rechazando un lanzamiento multiproceso. La entrada habitual `./grial caso.par`, incluido el build MPI de un proceso, conserva el camino previo y sus formatos; para evolucionar varios procesos se exige `--distributed` explícito. `--parallel-check` conserva su función de diagnóstico del runtime.

Para verificar configuración y reparto sin escribir resultados:

```bash
OMP_NUM_THREADS=2 mpirun -n 4 ./grial_mpi --distributed --check caso.par --grid 2 2
```

El arranque frío de `--check` no aloja campos. No es una reserva de recursos ni garantiza que una escritura futura sea posible.

## Parar y continuar por tiempo real

Se conservan `walltime_limit_seconds`, `walltime_margin_seconds` y `walltime_checkpoint_seconds` del `.par`. Las opciones `--walltime` y `--margin`, en segundos, sustituyen sólo el presupuesto de la sesión, sin cambiar los parámetros físicos ni el calendario persistente. Por ejemplo:

```bash
OMP_NUM_THREADS=2 mpirun -n 4 ./grial_mpi --distributed caso.par \
  --grid 2 2 --walltime 82800 --margin 900
```

El reloj incluye el arranque del programa/MPI/configuración. Cada rango comprueba su propio tiempo; una petición de cualquiera detiene todos en una frontera RK completa. Se descuentan el margen, el máximo costo de paso observado y el costo de checkpoint observado. Se guarda el estado físico y del observador antes de publicar `stopped_for_walltime`. La inicialización y la escritura también consumen tiempo: el margen debe calibrarse para la malla y sistema de archivos; no es una garantía frente a una caída abrupta o una primera escritura excepcionalmente lenta. Esta entrada no consulta automáticamente el deadline Slurm.

La petición de checkpoint por reloj se reduce entre procesos **antes** de entrar en las rutinas colectivas. Los relojes locales pueden cruzar el intervalo en pasos distintos; no se permite que sólo una parte de los rangos llame al escritor.

La continuación usa la generación exacta indicada por el recibo final, con un reparto nuevo si se desea:

```bash
OMP_NUM_THREADS=2 mpirun -n 8 ./grial_mpi --distributed --resume \
  data/mp5_hllc/mpi_checkpoints/step_000001234 \
  --grid 4 2 --walltime 82800 --margin 900
```

El directorio debe contener C01 y el cierre `OBSERVED` de C03; configuración, campos U/P, contadores, bandera de perturbación, presupuesto inicial/integrales, calendario e historias se validan y recuperan juntos. No se exige un `.par` nuevo ni se permiten cambios físicos implícitos al continuar. El tiempo objetivo se conserva desde el checkpoint. La ruta no importa aún checkpoints seriales V7/V8/V9 ni generaciones C01 sin observador; esas APIs de C02 siguen disponibles por separado.

```bash
OMP_NUM_THREADS=2 mpirun -n 8 ./grial_mpi --distributed --resume \
  data/mp5_hllc/mpi_checkpoints/step_000001234 --grid 4 2 --check
```

Este `--check` de recuperación sí aloja los campos y verifica su lectura completa e historias, pero no escribe resultados. Un checkpoint antiguo con historias adelantadas se rechaza: no se truncan ni revierten automáticamente. El formato de checkpoint es nativo y requiere un entorno binario compatible, como antes.

## Salidas y exclusión de escritores

Para `output_folder=data` y esquema `mp5_hllc`:

```text
data/mp5_hllc/
  run_status.dat
  mpi_global_diagnostics.dat
  mpi_conservation_budget.dat
  mpi_numerical_diagnostics.dat
  mpi_ppi_modes.dat                 # cuando se activa PPI
  mpi_GW_signal.dat                 # cuando se activa extracción
  mpi_m_dot.dat                     # cuando se activa extracción
  mpi_perturbation_event.dat        # si se aplica la perturbación
  mpi_checkpoints/step_000001234/    # C01 + observador C03
  mpi_fields/step_0000001234/       # piezas VTU + índice PVTU/PVD C04
```

`run_status.dat` es escrito sólo por raíz, mediante temporal, cierre, comprobación de tamaño y renombrado. Registra estado, generación, tiempo/paso/objetivo, reparto, paso inicial de sesión, duración y costos observados. `completed` requiere llegar al objetivo y publicar el checkpoint final; `stopped_for_walltime` es una parada limpia recuperable, no el fin científico. Un aborto puede dejar el último recibo en `running`; un exit0 aislado tampoco sustituye la aceptación numérica.

Un arranque frío rechaza un recibo existente y nunca reemplaza las generaciones de campo/checkpoint. El directorio exclusivo `.mpi_run_lock` evita dos escritores cooperantes sobre el mismo esquema. Se retira al cerrar correctamente; si persiste tras una interrupción, exige inspección manual del proceso y de la consistencia de la historia antes de retirarlo. No hay borrado/reintento automático. El lock no evita modificaciones ajenas de los archivos ni hace atómica toda la aplicación.

Las series se continúan desde el cursor C03, con épocas estrictamente crecientes y un solo evento de perturbación. La VTK se guarda en pasos físicos, sin duplicar el frame recuperado ni alterar dt para acertar tiempos nominales. Puede añadirse un frame terminal fuera del calendario al parar por reloj.

Para formar la colección temporal, con el escritor ya detenido:

```bash
python tools/parallel_vtk.py data/mp5_hllc/mpi_fields \
  --collection data/mp5_hllc/mpi_fields/frames.pvd
```

## Organización y errores

`parallel_application.f90` conecta las APIs existentes. `agree_command_line` acuerda argumentos antes de que los procesos elijan una ruta del programa. La lectura fría contrasta longitud/CRC64 del `.par` entre rangos antes de aplicar presets; no es un hash criptográfico ni admite edición concurrente del archivo. Los errores de lectura/validación en `parameters` y los de configuración/arquitectura en `initialization` llaman a `parallel_abort`, que aborta el comunicador activo. Fuera de una ejecución MPI inicializada termina con error, como las regresiones anteriores. No se añade MPI dentro de bucles OpenMP ni se modifica la recuperación, reconstrucción, flujo, CFL o RK.

## Verificación y alcance

```bash
make test-application
make USE_MPI=1 test-application-mpi
```

La suite corta entra desde `grial`, no desde un driver sustituto. La aceptación ampliada usa builds sin/conMPI opt/debug con bounds/traps; compara contra el driver C03/C04 y el ejecutable previo, continuaciones con otro grid y parada real con progreso, series/campos/VTK y rechazos de configuración/argumentos/generación. Los intentos se conservan en el expediente.

Este incremento integra la aplicación local. No demuestra escalabilidad, estabilidad física tardía de toros, miles de rangos o continuidad bajo Slurm/24h. Los lanzadores `run_mpi_segments.py` existentes usan el driver C03 y sus recibos; no se presentan como controladores automáticos de esta entrada nueva. La integración y calibración de relanzamiento en asignaciones reales pertenecen a V04. V01–V04 continúan pendientes hasta sus cierres específicos. No se alteran fuentes/campañas CADS congeladas ni la implementación magnética del padre.

La nueva ruta publica las series distribuidas descritas; no reproduce el texto de ETA de la ruta serial ni su archivo específico de normas analíticas de convergencia. Los balances D01 y los campos/checkpoints sí están disponibles para el análisis. Esas salidas seriales continúan disponibles por la entrada habitual.

## Resultado de aceptación

78 continuaciones con U/P exactos y70arranques fríos con error máximo 9.721390e-15; presupuesto≤3.784242e-14 y series≤1.033389e-14, criterio1e−11.78rechazos coordinados. Cuatro builds sin/conMPI opt/debug aprobaron37programas Fortran y tres grupos de integración cada uno. Se comprobaron40checkpointsV9 completos con1hilo y40oráculos del driver; ocho controles de checkpoint periódico y cuatro recuperaciones con perturbación ya aplicada sin repetirla. Las continuaciones de reloj avanzan antes de parar.

Los arranques fríos optimizados pueden cambiar últimos bits al evaluar funciones sobre bloques distintos; los debug fueron exactos en esta matriz. Los campos de las continuaciones son exactos. La comprobación estricta del formato legado usa un hilo, igual que los cierres anteriores; con dos hilos se observaron últimos bits distintos en una suma atómica cruda de reparación, sin cambio de campos ni contadores. No se cambió el umbral1e−11. Los intentos y correcciones de enlace/harness se conservan en DEVELOPMENT.md y verification_v1/v2/v3 del expediente.

## Ejecutable local disponible

`mpi/grial_mpi` apunta al binario optimizado MPI de la matriz aceptada, copiado a `.build/control_c_20261027_verified/grial_mpi`. El enlace anterior y su hash se conservan en `LOCAL_LINK_BEFORE.json` del expediente; su archivo `.build/mpi/grial_mpi` no se reemplazó. Para volver sólo al binario anterior, desde `mpi/`: `ln -sfn .build/mpi/grial_mpi grial_mpi`. La siguiente compilación habitual con `make USE_MPI=1` vuelve a gestionar el enlace mediante `.build/mpi`. La fuente previa está en `before/` y el cambio documentado en `cambios.patch`; no se revirtió nada automáticamente.
