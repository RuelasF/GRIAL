# C02: configuración, cambio de reparto e importación serial

**Cerrado localmente el 25/oct/2026.** [Evidencia verificada](../.validation/c02_restart_20261025_v1/VERIFIED.json).

Incremento del proyecto hidrodinámico CPU en `mpi/`. Conserva el formato de [C01](CHECKPOINT_DISTRIBUIDO_C01.md) y añade lectura de configuración antes de reservar la malla, recuperación con otro reparto radial/polar e importación de checkpoints seriales V7/V8/V9. La entrada multiproceso de producción mantiene su guarda. Los controles usan drivers separados.

Expediente: [fuentes, intentos y resultados](../.validation/c02_restart_20261025_v1/). El cierre se acredita mediante `VERIFIED.json`; las carpetas parciales conservan sus logs.

## Secuencia de uso

Las llamadas son colectivas, fuera de OpenMP y entre pasos RK completos. Todos deben acordar camino, operación, presencia de `status` y grid de procesos. Se necesita almacenamiento compartido. El runtime MPI debe estar inicializado.

```fortran
use parallel_checkpoint, only: read_parallel_checkpoint_configuration,load_parallel_checkpoint
! Sin U/UP/P reservados y con D01 inactivo:
call read_parallel_checkpoint_configuration(path,layout,status,process_shape=[nr,ntheta])
! Si status /= 0, detenerse o volver a cargar una configuración válida.
! Preparar coordenadas locales, métrica, fronteras, arrays de evolución y halos con layout.
call load_parallel_checkpoint(path,layout,status,redistribute=.true.)
! Sólo tras el éxito, iniciar explícitamente una nueva línea base D01.
```

El lector de configuración restaura dimensiones globales, parámetros físicos y numéricos, pisos efectivos, dato inicial FD/FM, perturbación, muestreo y controles de salida/tiempo real. Recalcula `g1` y las banderas de métrica. Construye el nuevo `domain_layout`; el reparto predeterminado es `[parallel_size,1]`, con phi completo. La configuración guardada debe ser válida según `validate_parameters`.

La fase de configuración aún no tiene campos físicos que reemplazar. Si devuelve error después de interpretar algunos registros, puede haber modificado parámetros globales: **no se debe seguir con esos parámetros parciales**. Los errores de lectura se coordinan mediante el argumento opcional de estado añadido al lector de metadatos compartido. Los errores semánticos del validador existente detienen la ejecución; no constituyen una configuración recuperable.

La fase de estado comprueba configuración, bloques, época y contenido en temporales locales antes del acuerdo colectivo que reemplaza `UP/U/P`, contadores y época. Ante un rechazo previo a ese acuerdo, preserva los campos, halos, contadores y época existentes. D01 debe estar finalizado al cargar. La llamada histórica `load_parallel_checkpoint(path,layout,status)` sigue exigiendo el mismo reparto.

No hay overrides generales de configuración en este incremento: la recuperación compara los bytes canónicos, incluidos los nombres de salida. El driver de pruebas abre directorios nuevos después de recuperar. La política de segmentos y sus rutas corresponde a C03.

## Lectura por intersecciones

El manifiesto se valida primero contra su propio número de rangos y su grid guardado. Cada descriptor debe corresponder a `build_domain`: forma, offset, bytes esperados y checksum. Después se compara la forma global, halo y periodicidad con el destino. No se exige igualdad entre el grid anterior y el nuevo cuando se solicita `redistribute=.true.`.

Cada archivo anterior se verifica completo una vez entre los rangos nuevos: el verificador es `old_rank mod new_size`. Cada rango lee las filas físicas de `UP/P` cuya intersección global le pertenece. Usa posiciones de archivo int64 y el orden Fortran `[x,y,z,variable]`. Se comprueba que cada celda local tenga cobertura única, junto con finitud, densidad/presión positivas y causalidad. Raíz lee también las cabeceras necesarias para sumar la historia de contadores en orden de rango anterior.

No se reúne la malla en raíz ni se reserva un bloque anterior completo. Los temporales de estado son dos arrays de tamaño local y una máscara bidimensional de cobertura. El manifiesto crece con el número de bloques anteriores. Se vuelve a generar la métrica estática y se rellenan halos desde las primitivas guardadas; la recuperación distribuida no introduce otra inversión C2P.

El checksum y los bytes de configuración mantienen su contrato C01. Se asume que una generación publicada es inmutable durante su lectura; no se ofrece recuperación de procesos MPI perdidos ni resistencia a una modificación concurrente del almacenamiento.

## Convención de contadores históricos

Con el mismo reparto se preservan los contadores y acumuladores locales originales. Si cambia el grid o el número de procesos, sus sumas históricas se asignan a rango 0 y los demás empiezan en cero. Se conserva el total global, con comprobación de desbordamiento int64 y FP64 antes de sumar. Esta convención se anuncia en la salida.

Ese total es historia agregada, no un mapa espacial ni una atribución al subdominio nuevo de rango 0. El trabajo posterior se acumula en su rango ejecutor. Los contadores de evaluaciones de caras pueden cambiar con el nuevo reparto por el trabajo de interfaces; los eventos físicos agregados deben coincidir. Los acumuladores crudos `numerical_delta_u` se preservan como suma, sin reinterpretarlos como un presupuesto integrado. La continuidad del presupuesto ponderado y sus máscaras queda en C03.

## Importación de V7/V8/V9

```fortran
use parallel_checkpoint, only: read_serial_checkpoint_configuration,import_serial_checkpoint
call read_serial_checkpoint_configuration(serial_path,layout,status,process_shape=[nr,ntheta])
! Comprobar status y preparar malla/métrica/fronteras/halos según layout.
call import_serial_checkpoint(serial_path,layout,status)
```

El primer lector indexa los registros secuenciales y sus cierres sin reservar el campo global. Raíz calcula tamaño y CRC64 del archivo con un buffer fijo; cada rango interpreta metadatos y lee después únicamente sus filas conservadas. Se comprueban nuevamente tamaño/checksum antes y después de preparar el estado. El archivo fuente no se modifica.

Estos formatos no guardan primitivas ni el último `dt`. Se recupera `P` mediante el mismo C2P del arranque serial, sobre una copia local por celda; `UP` se conserva y `U` se sincroniza con él, como en la ruta histórica. El trabajo de esta inversión inicial no incrementa los contadores numéricos. `dt=0` indica que no estaba guardado y el siguiente paso lo calcula de nuevo. Se restablecen tiempo, paso y bandera de perturbación.

| Origen | Información conservada o valor heredado |
| --- | --- |
| V7 | Metadatos físicos/numéricos e iniciales; pisos heredados `rho=1e-10`, `p=1e-13`; contadores y acumuladores cero. No se puede recuperar una historia que no estaba escrita. |
| V8 | Añade controles de ejecución, pisos efectivos y cierre verificable; sin contadores persistentes, que empiezan en cero. |
| V9 | Conserva además contadores y acumuladores crudos; el total histórico se asigna a raíz. |

Para V7 se fijan explícitamente los valores de arranque ausentes: `checkpoint_interval=100`, `extraction_stride=10`, `mass_monitor_stride=1000`, diagnósticos numéricos desactivados, VTK físico, límite de tiempo real cero, margen 900 s y estimación de checkpoint cero. El lector histórico anuncia la ausencia de pisos/contadores. La elección de arquitectura recompone el halo necesario para el esquema/sensor actual, siguiendo el arranque serial existente. No se inventan archivos anteriores ni un presupuesto D01 persistente.

La importación admite **registros nativos con marcadores de 4 bytes**. Rechaza registros segmentados, versiones distintas, truncados, cierres/dimensiones incoherentes y registros sobrantes. No establece portabilidad entre endianness/ABI diferentes. V7 carece de cierre/checksum históricos: se verifican estructura, configuración y estado, pero el CRC calculado al importar no demuestra integridad anterior a esa lectura. No es una migración de formato autenticada.

## Verificación reproducible

```bash
make test-restart-parallel
make USE_MPI=1 test-restart-mpi
```

`test_parallel_restart` y `test_restart_mpi.py` incluyen 34 variantes: onda, Sod, KHI, Jet cilíndrico, Michel y toros FD/FM a=0/0.9, sagitales/3D, WENO5–HLLE y MP5–HLLC con sensor, y tres perturbaciones. Dimensiones desiguales: onda/toro `37×19×8` (sagital con `nz=1`), Sod/Michel `73×1×1`, KHI/Jet `37×1×19`. Son evoluciones cortas de ocho pasos, con `dt≤0.002`; no prueban estabilidad física tardía.

El control produce referencias independientes con un rango, guarda antes y después de la perturbación y exige igualdad exacta de `UP/P`, tiempo/paso/dt y bandera tras continuar. Recorre `1→2→4→2`, `4→8→1` y grids radiales/polares `2×2`, `1×4`, `4×2`, `2×4`. Los casos con una sola celda polar usan reparto radial. Los contadores iniciales superan 2^32. Se comprueba historia agregada y residuo de un **presupuesto nuevo** menor que 1e-11.

Las importaciones de cada versión se comparan con la ruta serial anterior (`read_checkpoint_metadata`/`read_checkpoint_state` y C2P de arranque), no consigo mismas. Los fixtures V7/V8 se construyen eliminando de registros V9 las extensiones ausentes y adaptando sus marcadores; no son archivos recuperados de una campaña histórica. Cada versión se importa con 1/2/4/8 rangos y se continúa también desde el checkpoint distribuido resultante.

Los rechazos incluyen corrupción/configuración incompleta con checksums válidos, bloques ausentes/truncados, NaN, rango incorrecto, configuración diferente, desbordamiento de historia y metadatos semánticamente inválidos; también V7/V9 truncados, versión no admitida, registros sobrantes y fuente alterada después de cargar su configuración. Los errores de estado deben dejar la memoria previa intacta.

La aceptación completa exige ambas matrices optimizada/debug, regresión C01, cuatro suites sin/con MPI, controles del programa habitual frente a binarios pre-C02 con iguales flags y todos los registros V9 idénticos. El expediente conserva el primer fallo del fixture plano: tenía masa de agujero negro 1 en Minkowski; se corrigió a 0 sólo en el test nuevo. No se relajaron tolerancias ni se cambió un kernel físico.

El primer `make test` MPI optimizado encontró además una fluctuación del control histórico de tiempo real: el arranque tardó 0.363 s y agotó el límite de 0.30 s antes del primer paso. El programa produjo correctamente un checkpoint en paso cero, pero ese subcontrol exigía avance positivo. Se conserva `mpi_opt_suite.failed_attempt1.log`; el control aislado pasó con los mismos parámetros en `walltime_retry.log`. No se alteraron el test antiguo, sus límites ni el solver. Este límite subsegundo es sensible a la carga de la máquina.

## Resultado del cierre

Configuración antes de alojar la malla, recuperación por intersecciones con otro grid radial/polar e importación serial V7/V8/V9. 1564 comparaciones opt/debug, campos/época/bandera y continuación exactos; residuo de presupuesto nuevo≤4.941e-14 frente a1e−11. 40 rechazos C02 sin reemplazo parcial de campos; regresión C01 de544comparaciones,54corruptos,44usos inválidos y12fallos reales de escritura. 35Fortran+dos grupos+8autopruebas porbuild en cuatro builds;120controles habituales con registrosV9 completos idénticos. Historia agregada en raíz al cambiar reparto; configuración fallida puede dejar parámetros parciales y exige detener/releer. El primer fixture plano y un fallo por arranque que agotó el límite subsegundo del control antiguo quedan documentados; no se cambiaron tolerancias ni ese control. Referencia41 y193archivos protegidos intactos. 
La matriz completa y las suites aceptadas están en `verification_v1`; el primer fallo de suite conserva su log separado. Se usaron GNU Fortran14.2.0/OpenMPI5.0.7, FP64, `ARCH_FLAGS` vacío y dos hilos por rango; debug añade bounds/traps invalid/zero/overflow. Los controles nativos usan un hilo y binarios pre-C02 de iguales flags. Fuentes previas/después, parche, parámetros, binarios, hashes, recibos y controles aceptados permiten revisar o revertir el incremento.

## Alcance y siguiente tarea

C02 restaura estado físico y configuración entre repartos. **C03** debe conservar líneas base/máscaras/mapas, presupuesto integrado, estado de observadores y calendario de series, y conectar la parada colectiva por tiempo real y los segmentos MPI. **C04** conserva la visualización distribuida. No se ha habilitado la entrada MPI de producción, ni cerrado continuidad entre asignaciones, ejecución multinodo o escalabilidad.

Las pruebas locales con ocho rangos y dos hilos por rango sobresuscriben seis cores físicos: verifican corrección, sin interpretar rendimiento. No se consultan ni cambian jobs CADS durante este incremento, ni sus fuentes congeladas, el padre magnético o la referencia hidrodinámica.

Graphify reutiliza el grafo existente de 906 nodos; la consulta de configuración/reinicio localizó seis nodos históricos con presupuesto de 1500 tokens y sin nueva inferencia LLM. Los módulos MPI recientes se verifican directamente en los fuentes, porque ese grafo no los representa aún.
