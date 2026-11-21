# C03: observadores y segmentos MPI

**Cerrado localmente, 26/oct/2026.** El cierre se acredita con [VERIFIED.json](../.validation/c03_observer_20261026_v1/VERIFIED.json).

Incremento CPU/FP64 de la copia hidrodinámica `mpi/`. Conserva el formato físico C01 y la configuración/redistribución C02. Añade continuidad de integrales D01, calendarios y series D02/D03 mediante `parallel_observer.f90`, un driver separado y un lanzador local. No retira la guarda multiproceso de `grial`, no envía jobs y no modifica el padre magnético, referencia ni campañas congeladas.

## Estado conservado

`parallel_diagnostics:diagnostic_history` exporta las sumas iniciales, seis flujos de frontera, fuentes, reparaciones ponderadas RK y perturbación: 50 FP64 por rango para las cinco conservadas. Guarda también pasos acumulados y trabajo adicional en interfaces MPI. Las conservadas están densitizadas; los factores de integración mantienen el contrato D01.

Con el mismo grid se recuperan las sumas locales originales. Al cambiar grid, se suman las contribuciones históricas y se asignan a raíz, con cero en los otros rangos; las contribuciones nuevas se acumulan donde se calculan. Se conserva la línea base inicial. Las reducciones con otro reparto pueden introducir redondeo, por lo que se comprueban presupuestos con el criterio previo 1e-11, junto con identidad de campos cuando corresponde.

El contador de evaluaciones adicionales ya realizadas se conserva separado de las realizadas después del reinicio. Para un tramo con n pasos y grid fijo, el incremento es 3*n*(caras_de_trabajo - caras_únicas). Aplicar el reparto nuevo a todos los pasos históricos sería incorrecto. Se comprueba la suma tramo a tramo.

Las máscaras geométricas, propietarios de caras, volumen y métricas se reconstruyen desde el layout nuevo. El núcleo MPI anterior no tenía mapas celulares persistentes de reparaciones ni un campo de densidad inicial de referencia: C03 no los inventa ni afirma restaurarlos. PPI y proxy de ondas son integrales instantáneas; desenrollado de fase, ajustes y Fourier se realizan en análisis posterior de las series.

`observer_calendar` conserva próximo tiempo nominal de salida y checkpoint, índice nominal de salida y últimos pasos muestreados para D01 y observables. El cursor de `output` conserva raíz/rutas, último tiempo de cada serie y último paso D01. Se mantienen `diagnostic_stride` y `extraction_stride`, sin alterar dt para alcanzar muestras. Cada checkpoint/parada fuerza una muestra terminal si aún no existe; puede añadir filas fuera de la cadencia regular, sin duplicar tiempos al continuar. Se comparan las muestras comunes con la continua.

## Generación completa

Una carpeta C03 contiene los archivos C01 (`configuration.bin`, `manifest.bin`, un bloque UP/P por rango y `COMMITTED`) y dos archivos nuevos:

- `observer.bin`: ABI/forma/grid anterior/época, calendario, cursor, historia D01 por rango y tamaños/CRC64 de siete historias: tres D01, PPI, ondas, acreción y evento de perturbación.
- `OBSERVED`: marcador final que liga los bytes de `COMMITTED` y `observer.bin` mediante tamaño/CRC64, con cierre verificable.

Se publica primero C01. Raíz serializa el estado del observador en memoria, escribe/cierra, relee tamaño/CRC y los contrasta contra los bytes originales. Sólo después renombra `OBSERVED.tmp` a `OBSERVED`. Un `COMMITTED` físico sin `OBSERVED` no basta para declarar una generación C03 recuperable. La generación anterior se preserva.

El sidecar ocupa aproximadamente 2360+400*rangos bytes en el ABI actual; no reúne la malla. Guardar reúne sólo los 50 números por rango. El lector mantiene esa tabla de metadata en cada rango nuevo y lee los campos con C02 por intersecciones locales. Este acceso a metadata y sus checksums, los bloques por rango y el costo de E/S deben medirse antes de producción con miles de procesos. El formato es nativo del ABI/endianness comprobados, sin promesa de portabilidad entre compiladores diferentes.

## Secuencia de recuperación

Las llamadas son colectivas, fuera de OpenMP y entre pasos RK completos, en almacenamiento compartido. Deben coincidir camino, operación y presencia de `status`. Se requiere un runtime inicializado y D01 inactivo antes de cargar.

```fortran
use parallel_checkpoint, only: read_parallel_checkpoint_configuration
use parallel_observer, only: observer_calendar,load_observed_checkpoint
! Sin UP/P/U reservados:
call read_parallel_checkpoint_configuration(path,layout,status,process_shape=[nr,ntheta])
! Comprobar status. Preparar malla local, métrica, datos de frontera y halos.
call load_observed_checkpoint(path,layout,calendar,status)
! Comprobar status antes de evolucionar; D01 queda activo sólo tras éxito.
```

El lector C02 suministra configuración y época. C03 contrasta su cabecera con esa época y valida marcador, sidecar, ABI, dimensiones, calendario, finitud, cursor y archivos de historia antes de llamar al recuperador físico C02. También comprueba tiempos estrictamente crecientes, último tiempo/paso y evento único. Un rechazo de C03 previo a la carga física deja campos/halos, contadores, reloj, calendario, cursor y activación D01 intactos. C02 conserva su propia recuperación transaccional ante bloques/configuración inválidos.

Después del éxito se inicializa D01 con la historia guardada, recalculando geometría/propietarios para el layout nuevo; se restauran cursor y calendario. El recuperador no trunca ni repara historias automáticamente. Si los archivos ya avanzaron más allá del checkpoint, se rechaza: una recuperación de una generación anterior requiere preparar explícitamente otra copia de las historias coherente con ese punto. Se asume que generaciones/historias no cambian concurrentemente durante su lectura.

Las APIs C01/C02 anteriores mantienen su alcance. Un V7/V8/V9 o un C01 antiguo sin sidecar puede importar campos con C02 e iniciar una línea base nueva; no contiene una historia D01/C03 que pueda recuperarse retrospectivamente.

## Guardado y muestras

```fortran
call initialize_parallel_diagnostics(layout)
call initialize_observer_calendar(calendar)
call sample_parallel_observer(calendar,force=.true.)
! Evolución con CFL global, halos y RK distribuido, e inyección D03.
call sample_parallel_observer(calendar)
! En una frontera de paso, con el calendario actualizado:
call save_observed_checkpoint(path,layout,calendar,status)
```

El guardado requiere directorio de generación nuevo e historias coherentes. Con `status`, E/S o generación incompleta devuelve estado distinto de cero; sin él, aborta coordinadamente. Uso inválido (D01 activo al cargar, desacuerdo entre rangos, parcial RK) aborta coordinadamente. El guardado por sí solo no fuerza muestra: el driver lo hace y actualiza el calendario antes de publicar.

## Driver y lanzador local

Compilar desde `mpi/`:

```bash
make -j4 USE_MPI=1 ARCH_FLAGS='' mpi-segment-driver
```

El ejecutable queda en `.build/mpi/tools/mpi_segment_driver`. Se usa una entrada separada porque la integración de producción, errores/salida C04 y validación en asignaciones reales siguen pendientes. Para una prueba se necesita un `.par` hidro frío, `numerical_diagnostics=true` y strides positivos. El lanzador fija una carpeta `data` dentro de un directorio nuevo y preserva el archivo original. El caso elegido debe caber en el grid r×theta; phi permanece completo en cada bloque.

```bash
python3 tools/run_mpi_segments.py \
  --exe .build/mpi/tools/mpi_segment_driver --par /ruta/caso.par \
  --run-dir /ruta/nueva/prueba_c03 --grid 2 1 --threads 2 \
  --segment-seconds 300 --margin-seconds 30 --max-segments 3
```

Para continuar explícitamente, incluso con otro reparto compatible:

```bash
python3 tools/run_mpi_segments.py \
  --exe .build/mpi/tools/mpi_segment_driver --par /ruta/caso.par \
  --run-dir /ruta/nueva/prueba_c03 --grid 2 2 --threads 2 \
  --segment-seconds 300 --margin-seconds 30 --max-segments 3 --resume
```

`--launcher-arg` permite opciones concretas del mpirun instalado; en local se ha usado `--launcher-arg=--oversubscribe --launcher-arg=--bind-to=none` para comprobar corrección con ocho rangos. No interpreta rendimiento. `--serial` permite el build sin MPI sólo con grid1×1. `--test-segment-steps` es una parada determinista para regresiones; cero la desactiva. Su solicitud se origina sólo en el último rango y se acuerda con todos, probando propagación colectiva. Los argumentos no se ejecutan a través de un shell.

Cada sesión reinicia reloj real, estimaciones de duración del paso y checkpoint. Conserva época física, bandera de inyección, contadores, integrales y muestreo. Comprueba relojes locales y coordina cualquier solicitud de parada al terminar un RK, reservando margen y estimación de escritura. El margen debe calibrarse en el sistema de destino: no garantiza sobrevivir a un SIGKILL, caída de nodo, almacenamiento o un paso más lento que lo previsto. No es recuperación de procesos perdidos.

El presupuesto empieza al entrar al driver e incluye su inicialización MPI/malla; el tiempo del controlador o la compilación queda fuera. En Slurm se debe proporcionar el tiempo restante de la asignación. `--max-segments` encadena procesos locales con presupuestos por sesión; la continuidad entre asignaciones y su integración con los límites del scheduler corresponde a V04.

El lanzador preserva logs y recibos por sesión, valida los enlaces de metadata, calcula SHA256 por bloques de lectura y registra hashes de todos los archivos de la generación. Antes de continuar compara ejecutable, parámetros originales/efectivos y generación inmutable. El lector Fortran verifica además CRC de cada bloque e historias. No relanza sobre un fallo, recibo ambiguo o falta de progreso. Un lock impide dos controladores simultáneos sobre el mismo directorio; un lock dejado por interrupción exige inspección manual.

| Retorno del lanzador | Significado |
|---|---|
| 0 | Tiempo físico objetivo alcanzado y generación coherente; o cadena ya completada verificada. |
| 2 | Cuota de segmentos agotada tras una parada limpia con progreso y checkpoint válido; requiere `--resume`. |
| 1 | Fallo de sesión, archivos cambiados, uso inválido o recibo inconsistente; revisar logs, sin continuación implícita. |

Se publican `mpi_segments.json`, `session_000001.log/.dat`, etc. El estado `running` no significa completado. Los checkpoints se guardan en `checkpoints/step_NNNNNNNNN`; no se sobreescriben. El formato de visualización de campos distribuidos queda para C04.

## Verificación

```bash
make test-observer test-mpi-segments
make USE_MPI=1 test-observer-mpi
```

La matriz usa ondas/Sod/KHI/Jet/Michel y FD/FM 3D/sagitales, espines0/0.9, WENO5–HLLE/MP5–HLLC y tres perturbaciones. Incluye repartos desiguales de1/2/4/8 procesos; donde ny=1 sólo se reparte radialmente. Se comparan todos los UP/P con lectura Python independiente y presupuestos/muestras comunes con 1e-11. La evolución se detiene/reanuda antes y después de la inyección; se exige evento único, línea base persistente y trabajo de interfaces sumado por tramo.

Los toros de continuidad fijan explícitamente rin=8M y rp,max=12M, como los controles de referencia, en Kerr–Schild: los presets de otras configuraciones no garantizan un toro cerrado al cambiar el espín. Estas mallas pequeñas y tiempos cortos comprueban persistencia/reparto, sin evaluar estabilidad física tardía.

La comparación PPI usa C0, Re/Im de Cm normalizados con C0, amplitud y cuerda de fase ponderada por amplitud: 2*min(Am1,Am2)*abs(sin(m*delta_phi/2)). Cada fase se comprueba además contra sus propios coeficientes. Esto evita exigir una fase absoluta precisa a un coeficiente casi nulo; un intento previo comparaba ángulos directamente con Am≈1e-11 pese a campos idénticos y se conserva como criterio defectuoso. El umbral de error complejo sigue siendo 1e-11; no se modifica el solver. Para las otras filas se registra la norma max(1,abs(fila)); los oráculos de las fórmulas de observables pertenecen a D02. El cierre C03 no calibra amplitudes de ondas ni certifica fases físicamente resueltas en todos los modos.

Se comprueban también parada real por reloj, cadena launcher1→2→4→2, cuota/completo, intento de reanudar campos modificados, controles del lanzador y recuperación rechazada sin mutación ante corrupción/falta/truncado de sidecar, marcador, campos o historias. Algunos sidecars corruptos tienen checksums recalculados para probar semántica. Un límite real de archivo se aplica después de inicializar MPI, sólo en raíz: debe impedir publicar el observador y dejar recuperable la generación anterior. El intento que aplicó ese límite demasiado pronto afectó a PMIx; queda conservado y no cuenta como prueba de E/S del solver.

El expediente contiene fuentes previas, parches, intentos y suites optimizadas/debug con bounds/traps, sin/con MPI. 448 comparaciones segmentadas/continuas, campos exactos; residuo persistente ≤3.262e−15 y error de series ≤5.499e−15 frente a 1e−11. 76 corrupciones rechazadas, cuatro fallos reales de escritura y cuatro cadenas con parada por reloj real; inyección única antes/después y cambio de reparto comprobados. 36 programas Fortran y dos grupos de integración en cuatro builds sin/con MPI optimizados/debug; 120 controles habituales con registros V9 completos idénticos. Padre/paquetes protegidos (193 archivos) y referencia (41) intactos. Estos son controles de continuidad, no validación tardía de PPI, escalabilidad, E/S de miles de rangos ni continuidad entre asignaciones reales. **C04, Control C de integración y V01–V04 continúan pendientes.**


## Resultados de aceptación — 26/oct/2026

[Recibo global](../.validation/c03_observer_20261026_v1/VERIFIED.json) y [matriz final](../.validation/c03_observer_20261026_v1/verification_v4/BUILDS_VERIFIED.json).

| Control | Resultado |
|---|---|
| Continuaciones de la matriz principal | 424: 184 por build MPI y 28 por build sin MPI |
| Cruce de la inyección y checkpoints periódicos | 24 comparaciones adicionales; evento único y calendario persistente |
| Error máximo absoluto de UP/P | 0 |
| Residuo máximo del presupuesto persistente | 3.2612e−15 |
| Error máximo normalizado de series | 5.4985e−15, criterio 1e−11 |
| Corrupciones de generación/historias | 76 rechazos |
| Límite real de archivo, parada por reloj y launcher | Cuatro controles de cada clase, uno por build |
| Regresión C02 | 1564 comparaciones y 40 rechazos; campos exactos |
| Regresión D01 | Matrices completas opt/debug aprobadas |
| Compatibilidad habitual frente a pre-C03 | 120 controles; campos, contadores, balances y registros V9 completos idénticos |

Los cuatro builds aprobaron `make test`: 36 programas Fortran, dos grupos de integración, ocho autopruebas del controlador de campaña, once del launcher MPI y cinco del launcher anterior. La cadena local 1→2→4→2 conserva la línea base de presupuestos; los controles adicionales cruzan la inyección con 1→4→8 procesos.

Se conservan los intentos fallidos: fixtures/aliases iniciales, fase casi nula mal comparada, límite de archivo aplicado antes de arrancar MPI, expectativa antigua del mensaje de rechazo D01 y un control heredado cuyo presupuesto subsegundo se agotó durante el arranque. Se corrigieron fixtures/harness o se repitió sin cambios el control temporal; no se cambiaron tolerancias físicas ni kernels. La evidencia definitiva está en `verification_v4`, `injection_boundaries` y `NATIVE_VERIFIED.json`.

Siguiente paso: C04, visualización verificable de bloques. Control C, entrada MPI de producción, continuidad entre asignaciones reales y escalabilidad permanecen pendientes. No se trasladó C03 a fuentes remotas congeladas.
