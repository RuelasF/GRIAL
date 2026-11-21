# C01: generación de checkpoint distribuido

**Cerrado localmente el 25/oct/2026.** Proyecto hidrodinámico CPU en `mpi/`: guardado y recuperación **con el mismo reparto y el contexto ya preparado**. [Resultados verificados](../.validation/c01_checkpoint_20261025_v1/VERIFIED.json), con fuentes previas, parche, parámetros, binarios, logs y hashes.

## API y punto seguro

```fortran
use parallel_checkpoint, only: save_parallel_checkpoint,load_parallel_checkpoint
call save_parallel_checkpoint(generation_path,layout,status)
! Tras preparar la misma configuración, malla, métrica y halos:
call load_parallel_checkpoint(generation_path,layout,status)
```

Las llamadas son colectivas y ocurren fuera de regiones OpenMP y entre pasos RK completos. El mismo camino, operación y presencia del argumento `status` deben acordarse entre rangos. La guarda de producción multiproceso de `grial` permanece activa; los controles usan drivers separados.

`status=0` significa éxito y `status=1` rechazo común de E/S, configuración incompatible o contenido inválido. Sin ese argumento, un error provoca rechazo colectivo. Un uso inválido del API —dominio incoherente, época no finita, llamada a mitad de RK— se rechaza colectivamente antes de escribir. No se transforma un fallo en una continuación automática.

El guardado permite D01 activo y comprueba que no haya una etapa pendiente. La recuperación exige **D01 finalizado** para evitar continuar con un presupuesto anterior incompatible. El llamador puede iniciar explícitamente un presupuesto nuevo después de cargar; eso establece otra línea base. Restaurar la historia y la línea base originales corresponde a C03.

## Archivos y publicación

Cada llamada de guardado utiliza un directorio nuevo, cuyo padre ya existe:

```text
generation_path/
  configuration.bin
  block_000000000.bin
  block_000000001.bin
  ...
  manifest.bin
  COMMITTED
```

Hay `n_ranks+3` archivos por generación. El directorio se crea exclusivamente; si ya existe, la operación falla sin sobrescribirlo. No se elige una generación reciente por su nombre o fecha: el lector necesita el marcador y valida todo el contenido. Los directorios incompletos se conservan como evidencia y se requiere otro nombre para repetir el guardado.

1. Cada rango genera los bytes canónicos de la configuración global y acuerda su tamaño/checksum, época y bandera de perturbación.
2. Raíz crea el directorio y escribe `configuration.bin`. Todos comprueban su visibilidad y contenido, por lo que se necesita un sistema de archivos compartido.
3. Cada rango escribe y cierra su bloque. Comprueba el tamaño esperado, calcula su checksum y vuelve a leerlo, contrastando campos y contadores con la memoria original. El acuerdo requiere éxito de todos los rangos.
4. Raíz reúne sólo los descriptores, escribe/cierra el manifiesto y comprueba su tamaño/checksum.
5. Raíz escribe `COMMITTED.tmp`, lo vuelve a leer y publica `COMMITTED` mediante un renombrado dentro del directorio. Éste es el punto de publicación de la generación.

Un bloque ausente o incompleto impide llegar al paso 5. Las generaciones anteriores se conservan sin modificaciones. La sincronización comprueba errores de procesos/E/S; no implementa recuperación de un proceso MPI perdido. Se asume el comportamiento POSIX del sistema compartido y no se ha demostrado durabilidad frente a corte eléctrico o caída completa del almacenamiento.

## Datos conservados y formato

El formato usa streams binarios nativos con versión, tamaños ABI, marcador de orden de bytes y cierre. No es un contenedor HDF5 ni un archivo V9 serial completo. Se comprueban FP64, enteros, formas y coordenadas de los bloques antes de reservar/cargar el estado. La portabilidad entre compiladores/representaciones diferentes no se da por demostrada.

`configuration.bin` reutiliza `output:write_checkpoint_metadata`, extraído del escritor V9 sin cambiar su orden de registros predeterminado. En C01 escribe dimensiones **globales** y una representación canónica: bandera de perturbación aplicada, contadores y acumuladores locales se ponen a cero sólo en esos bytes. Sus valores efectivos se guardan por separado. Incluye límites, métrica/espín/EoS, esquema/solver/sensor, parámetros iniciales FD/FM y de otros casos, pisos efectivos, perturbación, muestreo, salida y controles de tiempo real.

El lector C01 compara esos bytes con la configuración preparada por el llamador, incluidos los nombres/rutas de salida. No reconstruye ni cambia automáticamente la configuración. C02 añadirá la carga de configuración y los overrides permitidos, importación serial y otro reparto.

Cada bloque contiene:

- Versión/ABI, rango, número de procesos, dimensiones globales/locales, offsets, grid de procesos, halos y periodicidad.
- Tiempo físico, último `dt`, paso y bandera `perturbation_applied`.
- Contadores locales int64 y acumulador crudo local `numerical_delta_u`.
- Las cinco conservadas físicas `UP` y las cinco primitivas físicas `P`, en orden `[x,y,z,variable]`, y el cierre del archivo.

Se guardan las primitivas para no introducir una nueva inversión al recuperar, incluida una parada inmediatamente después de la perturbación. La métrica estática, coordenadas y datos de frontera se preparan desde la configuración; los halos se regeneran. La copia `U` se sincroniza con `UP`, pues el siguiente RK la vuelve a formar. RHS/fuentes y la capa extra no usada de índice `-nghost` no forman parte del estado persistente.

El manifiesto registra configuración, época y un descriptor por rango: rango, tres tamaños, tres offsets, tamaño de archivo y checksum. La tabla se contrasta con `build_domain`; no basta con que existan los nombres de archivo. Los checksums CRC64 se calculan con el polinomio `0x42F0E1EBA9EA3693`, valor inicial cero; el vector `123456789` produce `0x6C40DF5F0B497347`. Se usan para detectar corrupción accidental, junto con versiones, tamaños, cierres y validación física.

## Recuperación sin modificación parcial

Todos verifican el commit, el manifiesto, la configuración y su propio bloque. El lector comprueba ABI, reparto, offsets, época, banderas, contadores, finitud, densidad/presión positivas y causalidad. Reserva dos arreglos temporales locales y sólo después del acuerdo colectivo reemplaza `UP/P`, sincroniza `U`, restaura época/bandera/contadores y llena halos. Si un rango encuentra un error antes de ese punto, ninguno cambia el estado físico, los contadores o la época.

No se reúne la malla en raíz. Los temporales de campos crecen con las celdas **locales**; el manifiesto crece con el número de rangos. El tamaño del payload es aproximadamente `80*N_local` bytes por bloque, más cabecera/contadores. La lectura de comprobación añade E/S y memoria temporal al guardado. Su coste y la presión sobre metadatos con miles de archivos deben medirse en V02/V03 antes de elegir entre este formato y MPI-IO para producción.

## Verificación

```bash
make test-checkpoint-parallel
make USE_MPI=1 test-checkpoint-mpi
```

El programa `test_parallel_checkpoint` usa bloques desiguales y comprueba dos generaciones, rechazo de sobrescritura, recuperación al mismo reparto y continuación exacta. Incluye ondas, Sod, KHI, Jet, Michel y FD/FM con a=0/0.9 en 3D/sagital, WENO5–HLLE/MP5–HLLC, tres perturbaciones y contadores mayores que 2^32. Las evoluciones son controles cortos de E/S/reparto, no validación física tardía.

El harness Python inspecciona todos los archivos con una implementación independiente del checksum. Prueba falta de commit/bloque, corrupción de configuración/manifiesto/datos, truncado y contenido inválido con checksums actualizados. Los fallos de un rango se propagan sin cambiar la memoria de los otros. Provoca además un límite real de archivo de 4096 bytes: una variante devuelve el error de E/S y otra interrumpe un proceso. Ambas deben dejar la generación sin publicar y la anterior recuperable.

**Hallazgo de la prueba de E/S:** al ignorar SIGXFSZ, el runtime GNU Fortran pudo devolver éxito aunque dejara un archivo truncado. Por ello C01 no confía sólo en `iostat/close`: exige tamaño esperado y lectura de verificación antes del commit. El intento que detectó el problema permanece en `verification_v3`; el intento anterior de la señal también se conserva. Las correcciones del harness no relajan los umbrales físicos.

La aceptación exige suites completas sin/con MPI optimizadas/debug con bounds/traps, matriz de 1/2/4/8 rangos, regresión D03 y comparación habitual contra binarios pre-C01 con iguales flags. El formato V9 debe conservar todos sus registros por bytes. Fuentes previas, parches, logs y hashes quedan en el expediente.

## Resultado del cierre

La ejecución aceptada está en `verification_v5`. Las ocho particiones `1×1`, `2×1`, `1×2`, `4×1`, `1×4`, `2×2`, `4×2` y `2×4`, optimizado/debug, completaron **544 comparaciones** y verificaron 1088 generaciones. Los campos, primitivas y continuación física coincidieron exactamente. El residuo máximo del presupuesto nuevo fue **4.1733e−14**, frente al criterio 1e−11.

Pasaron **54 pruebas de corrupción**, incluidos NaN y rango incorrecto con checksums actualizados; **44 rechazos de uso inválido** y **12 fallos reales de escritura** (seis errores devueltos y seis caídas de proceso). Las generaciones fallidas quedaron sin commit y las anteriores siguieron recuperables. El lector definitivo rechazó además, en optimizado/debug y sin modificar memoria, la generación truncada de `verification_v3` que tenía commit y checksums coherentes con sus bytes incompletos.

Cuatro suites sin/con MPI, optimizado/debug, aprobaron **34 programas Fortran, dos grupos de integración y ocho autopruebas de campaña por build**. La regresión D03 repitió 864 comparaciones y 44 rechazos. En **120 controles del programa habitual**, todas las conservadas, épocas, contadores, balances, series y registros completos V9 fueron idénticos a los binarios pre-C01 de iguales flags. Se usaron GNU Fortran14.2.0/OpenMPI5.0.7 y seis cores físicos locales; ocho rangos con dos hilos son una prueba de corrección con sobresuscripción, sin interpretación de rendimiento.

Los 193 archivos protegidos, incluidos fuentes/binario del padre y ambos paquetes CADS, permanecieron intactos; la referencia completa de41archivos también se verificó. No se consultaron ni modificaron jobs remotos durante C01. Graphify reutilizó el grafo de906nodos para localizar las rutinas seriales históricas, con presupuesto de consulta1500tokens y sin nueva inferencia LLM; los incrementos MPI se revisaron directamente en los fuentes.

## Alcance posterior

C01 no habilita `restart=true` multiproceso, cambio de reparto, continuación de series, presupuesto D01 persistente, parada segmentada MPI ni visualización distribuida. Siguen **C02**, C03, C04 y la validación real en cluster. Las fuentes CADS, el proyecto magnético y la referencia inmutable permanecen separados.
