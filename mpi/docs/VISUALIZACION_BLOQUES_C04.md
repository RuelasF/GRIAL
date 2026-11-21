# C04: visualización de bloques distribuidos

**Cerrado localmente, 26/oct/2026.** El cierre se acredita con [VERIFIED.json](../.validation/c04_visualization_20261026_v1/VERIFIED.json). Esta salida se usa en el driver CPU experimental; la entrada MPI de producción conserva su guarda.

## Formato y contrato

`parallel_visualization.f90` ofrece una operación colectiva:

```fortran
call save_parallel_vtk(generation_path, layout, status)
```

Todos los rangos deben llamarla fuera de OpenMP y al terminar un paso RK, con el mismo directorio, configuración y época. Comprueba el layout contra la malla alojada, estado finito, dimensiones/offsets y acuerdo de tiempo/paso/geometría. Si D01 está activo, también exige que no haya un RK parcial. Los usos inválidos abortan coordinadamente; `status` distinto de cero señala un fallo colectivo de E/S. La rutina no modifica campos, contadores, tiempo o calendario.

El directorio padre debe existir. El escritor crea un directorio exclusivo por instantánea y rechaza uno ya existente. Cada proceso escribe sólo sus celdas físicas, sin ghosts; no hay Gather de puntos o campos ni un volumen global en raíz. La memoria temporal del escritor escala con el bloque local (puntos, conectividad y un escalar de trabajo), no con toda la simulación. Raíz escribe únicamente los índices pequeños de piezas y tiempo, de tamaño O(rangos).

Una instantánea completa contiene:

```
step_0000000123/
  piece_0000000000.vtu
  piece_0000000001.vtu
  ...
  frame.pvd
  frame.pvtu
```

Son archivos XML VTK auténticos: `.vtu` con geometría/celdas y datos binarios appended raw, `.pvtu` paralelo con referencias relativas a esas piezas, y `.pvd` con el tiempo físico de la instantánea. No se referencia VTK legacy desde un índice XML paralelo. [Especificación primaria VTK](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html).

Los escalares y coordenadas son Float64; conectividad, offsets e IDs son Int64, tipos de celda UInt8 y longitudes de payload UInt64, con orden de bytes BigEndian explícito. Raw appended contiene bytes binarios que un parser XML genérico no puede leer como un documento XML completo: se debe usar VTK o el lector C04, que separa cabecera y payload. No hay compresión ni dependencia VTK/HDF5 en el ejecutable Fortran.

Cada rango escribe a `.tmp`, cierra y contrasta el tamaño real con la posición final esperada antes de renombrar su pieza. Raíz publica `frame.pvtu` al final, sólo tras acuerdo de que todas las piezas y el índice temporal se cerraron correctamente. Un directorio sin ese índice no es una instantánea completa. Se conservan los archivos parciales para diagnóstico y la instantánea anterior permanece intacta; no se limpia ni reintenta automáticamente sobre un directorio fallido. No se promete atomicidad conjunta de VTK y checkpoint ni recuperación automática frente a fallos del sistema de archivos.

## Celdas, geometría y variables

- Se utilizan las caras numéricas verdaderas. Los ejes transversales colapsados se representan en su centro, igual que en la salida serial; x/r conserva sus caras incluso con una celda. Hay líneas 1D, cuadriláteros 2D y hexaedros 3D.
- En geometría esférica/esferoidal 3D, una cara sobre theta=0/pi se colapsa a un punto por radio y su celda se representa con un prisma triangular (`VTK_WEDGE`). Esto depende de la posición global del polo, no de que un bloque contenga todo theta. Un bloque interior no crea polos artificiales.
- Para un dominio azimutal completo, el último anillo se conecta al primero dentro de cada pieza polar. En un dominio azimutal parcial no se cierra esa frontera. Los puntos de interfaz radial/polar pueden repetirse entre piezas; las celdas globales pertenecen a un único bloque. No hay celdas fantasma que deban eliminarse al sumar datos.
- Se reutiliza `output::vtk_polar_position`: mapeo Kerr–Schild físico o visualización oblata sin giro radial, según `vtk_mapping_id`; radio logarítmico convertido a físico. Los polos se colapsan exactamente. La selección del mapeo sólo afecta a la figura.

Los siete campos celulares son los mismos de la salida serial: `Density`, `Pressure`, `Velocity_X/Y/Z`, `Lapse_alpha` y `Shift_beta_r`. Los nombres de velocidad son históricos: en geometría curvada representan componentes en la base coordenada correspondiente, **no un vector cartesiano transformado**. En radio logarítmico se convierten la velocidad y shift radiales multiplicando por r, igual que en VTK serial. No se añade una segunda recuperación de primitivas ni se cambia el estado nativo.

`Density` es la densidad propia primitiva rho, no D densitizada. El volumen geométrico de las celdas VTK sirve para visualizar su topología; no sustituye los pesos métricos del solver ni basta para reconstruir la masa relativista integrando rho. Para balances usar D01 o los campos conservados del checkpoint.

`GlobalCellId` es el índice global de celda, base cero y con x/r variando más rápido: `(i-1)+Nx*((j-1)+Ny*(k-1))`. `RankId` identifica al escritor. La metadata por pieza registra `TimeValue`, `Step`, `GlobalShape`, `LocalShape`, `GlobalOffset`, `ProcessGrid`, `Rank`, `MappingId`, `LogRadius` y `Spin`. Sirve para reconstruir la cobertura incluso con bloques desiguales.

## Calendario y reinicios

El driver `tools/mpi_segment_driver.f90` escribe una instantánea inicial, las salidas al cruzar el calendario `next_save` C03 y una final o terminal por parada temporal. El nombre usa el paso físico, no un contador local de sesión ni el número de procesos. No se modifica dt para alcanzar un instante nominal de salida: `TimeValue` y `.pvd` contienen el tiempo realmente alcanzado. Un paso ya escrito no se repite al coincidir salida periódica y final.

En una recuperación C04 se omite el estado de checkpoint ya publicado por la sesión anterior, y el calendario se restaura desde C03. Al cambiar el reparto, los nuevos frames pueden tener otro número de piezas; no cambian los IDs ni la forma global. La parada puede añadir una instantánea terminal fuera de la cadencia regular; las series continuas/segmentadas se comparan en sus pasos comunes.

El sidecar/formato C03 permanece sin cambios. Importar un checkpoint anterior a C04 permite continuar, pero no inventa las figuras históricas que nunca se guardaron. Si el proceso falla entre salida y checkpoint, no se retrocede ni sobrescribe una historia/instantánea adelantada: se debe preparar explícitamente una copia coherente con el punto de recuperación. El launcher actual sigue comprobando la identidad del ejecutable y parámetros al reanudar.

## Uso

```bash
cd ~/Codigo/GRIAL/mpi
make -j2 USE_MPI=1 ARCH_FLAGS='' mpi-segment-driver
```

Los argumentos del driver/launcher siguen documentados en [C03](OBSERVADORES_SEGMENTOS_C03.md). Las figuras se escriben dentro de `<output_folder>/<scheme_name>/mpi_fields/step_XXXXXXXXXX/`.

Para abrir una instantánea en ParaView se puede seleccionar su `frame.pvtu` (piezas conjuntas) o `frame.pvd` (con tiempo explícito). Para reunir la secuencia temporal, después de verificar que el escritor ya cerró esas instantáneas:

```bash
python tools/parallel_vtk.py data/weno5_hlle/mpi_fields/step_0000000123
python tools/parallel_vtk.py data/weno5_hlle/mpi_fields \
  --collection data/weno5_hlle/mpi_fields/frames.pvd
```

Abrir `frames.pvd` permite recorrer el tiempo. El segundo comando vuelve a generar atómicamente sólo ese índice de colección. Rechaza directorios `step_*` incompletos, tiempos/pasos duplicados o no crecientes y cambios de forma global, en lugar de ocultarlos. Los archivos de pieza y frame nunca se reescriben. Se asume que esos archivos no cambian concurrentemente durante la inspección.

`tools/parallel_vtk.py` usa sólo la biblioteca estándar. Lee una pieza cada vez y comprueba formato binario, longitudes/offsets, finitud, conectividad/tipos, IDs/propietarios, forma/reparto, cobertura y tiempo; contrasta el índice con las piezas y su `.pvd`. La inspección normal no aloja un campo global. Para un análisis que lo necesite explícitamente:

```python
import sys
sys.path.insert(0, 'tools')
from parallel_vtk import read_frame
frame = read_frame('data/weno5_hlle/mpi_fields/step_0000000123', reconstruct=True)
rho = frame['fields']['Density']  # array('d'), x/r varía más rápido
nx, ny, nz = frame['shape']
```

`reconstruct=True` aloja siete escalares globales en el proceso de análisis. Esto es una operación posterior optativa, no la estrategia de E/S del solver. El lector valida la estructura, no lleva checksums científicos: una alteración que conserve valores finitos y estructura puede no detectarse. Los checkpoints C01–C03 conservan su propio contrato de integridad y siguen siendo la fuente de reinicio; VTK no los sustituye.

## Verificación y límites

```bash
make test-visualization
make USE_MPI=1 test-visualization-mpi
```

El primer target también relee la salida con el inspector independiente. La matriz completa compara geometría/celdas por ID y campos reconstruidos para repartos desiguales de 1/2/4/8 procesos; ejes x, x–z, x–y, 3D, cilindro, EF y Kerr–Schild, r/log(r), polos completos y dominios angulares parciales, cortes sagital/ecuatoriales y ambos mapeos. El consumidor oficial VTK se usa en la aceptación local desde un entorno temporal aislado; no es dependencia obligatoria de la suite habitual. Se comprueba además la compatibilidad con VTK legacy usando sus siete escalares y los vértices de cada celda.

La conversión de v_radial/shift desde log(r) puede seleccionar kernels exp escalares/vectoriales diferentes según bucle/compilación. Su error normalizado se registra con el criterio numérico existente1e−11; los escalares sin esa conversión permanecen sujetos a identidad. Los campos nativos UP/P se comparan por separado con checkpoints, sin tolerancia adicional. El primer intento que exigía identidad de las conversiones queda conservado y medido en el expediente.

También se verifican colisiones, fallos de un escritor, límite real de tamaño de archivo aplicado después de arrancar MPI, argumentos inconsistentes y corrupciones estructurales; cadenas continuas/segmentadas SRHD, Michel y FD/FM con cambio de reparto. Las suites optimizadas/debug usan bounds/traps y conservan los intentos fallidos.

Esto no certifica rendimiento de E/S con miles de archivos/rangos, VisIt o una sesión GUI de ParaView, estabilidad tardía de toros, multinodo/GPU ni continuidad bajo Slurm. C04 no retira la guarda del programa principal. Después siguen Control C/integración de entrada y errores, V01–V04 y medición de escalabilidad.


## Resultado final — 26/oct/2026

234 controles de geometría/salida con el lector oficial VTK y 40 cadenas continuas/segmentadas; campos nativos exactos. Error de conversión frente a legacy ≈1.11e-16; presupuesto ≈1.44e-15 y series ≈2.208e-15 frente a 1e−11. 72 corrupciones de pieza/índice y ocho de colección rechazadas; diez usos colectivos inválidos, doce controles de E/S y cuatro límites reales de archivo. Cuatro builds sin/con MPI opt/debug: 37 programas Fortran, dos grupos de integración y autopruebas; C03 repetido con cadenas/tiempo real/corrupción y 120 controles habituales con registros V9 idénticos.

Los 234 controles de topología incluyen 60 comparaciones con los campos y vértices del VTK legacy. La geometría entre repartos es idéntica; longitudes/áreas/volúmenes de todas las celdas son positivos en VTK9.7.1. Los campos de densidad/presión y componentes sin conversión son idénticos; las diferencias registradas frente a legacy proceden sólo de convertir v_r/shift desde log(r).

Los intentos de desarrollo y la primera matriz se conservan. Se corrigieron fixtures y comparadores: lectura de todos los escalares legacy, redondeo de exp y reconocimiento de prismas equivalentes a hexaedros degenerados legacy en phi parcial. Un control R04 anterior agotó su presupuesto subsegundo durante el bootstrap singleton. Se reprodujo con los ejecutables anterior/actual; el harness MPI usa ahora mpirun de un proceso, como el arranque normal, sin modificar sus límites. Los kernels físicos y las tolerancias científicas no se cambiaron.

[Matriz final](../.validation/c04_visualization_20261026_v1/verification_v2/BUILDS_VERIFIED.json), [controles de colección](../.validation/c04_visualization_20261026_v1/collection_controls/VERIFIED.json), [compatibilidad nativa](../.validation/c04_visualization_20261026_v1/NATIVE_VERIFIED.json). Fuentes previas, hashes y parche en el expediente.193archivos protegidos y referencia41 intactos.

[Ejemplo de densidad desde los bloques](../.validation/c04_visualization_20261026_v1/comparacion_salida_MPI.png): FD/FM en una malla37×19×8, t≈0.065M, misma densidad entre continua de un proceso y reanudada2×2. Es una demostración corta de salida/reinicio; no corresponde a los pilotos3D prolongados de CADS ni valida estabilidad/PPI.

Siguiente: Control C, revisión conjunta de entrada/configuración/errores y conexión de estas rutas al programa de producción. La guarda permanece hasta verificar esa integración; V01–V04 y escalabilidad pendientes.
