# R04 — Continuación por segmentos y continuidad de diagnósticos

23/oct/2026 local, 24/oct UTC. Implementación y controles locales en `mpi/`; el solver sigue siendo OpenMP de un proceso. El punto de entrada del proyecto y el siguiente paso están en `CONTEXTO_MPI.md` y `tasks/todo.md`.

## Ejecutar y continuar

El lanzador usa Python3 con biblioteca estándar y llama únicamente al ejecutable local. No solicita ni envía trabajos Slurm. Se ejecuta desde `~/Codigo/GRIAL/mpi/`:

```bash
python3 tools/run_segments.py benchmarks/restart/michel.par \
  --run-dir .validation/michel_segmentos_nuevo \
  --segment-seconds 75 --margin-seconds 5 --checkpoint-seconds 25 \
  --threads 4 --max-segments 3
```

Para KHI usar `benchmarks/restart/khi.par`. El directorio debe ser nuevo: el lanzador evita sobrescribir otra ejecución. Copia el archivo original, conserva su SHA256 y el del ejecutable, y fija `output_folder=data` dentro del directorio elegido. Las opciones de tiempo real sustituyen únicamente los controles correspondientes del `.par`; no cambian la malla, reconstrucción, CFL ni tiempo físico final. Usa `OMP_NUM_THREADS` y desactiva ajuste dinámico de hilos.

Cada proceso atiende su margen antes de iniciar otro RK, termina con un checkpoint completo y publica su estado. El lanzador comprueba código de salida, estado terminal, avance físico, ruta del checkpoint y concordancia de paso/tiempo mediante `grial --check` antes de registrar la continuación. Valida de nuevo el checkpoint al empezar el siguiente segmento.

Códigos del lanzador:

- `0`: objetivo físico completado. No inicia más segmentos.
- `2`: se agotó `--max-segments`, con una parada planificada y checkpoint disponible.
- `1`: fallo, incoherencia o ausencia de avance. No continúa automáticamente.

Después de código 2, repetir el comando con `--resume`:

```bash
python3 tools/run_segments.py benchmarks/restart/michel.par \
  --run-dir .validation/michel_segmentos_nuevo --resume \
  --segment-seconds 75 --margin-seconds 5 --checkpoint-seconds 25 \
  --threads 4 --max-segments 3
```

`--max-segments` limita esta invocación. Puede ajustarse el presupuesto de un nuevo segmento. El archivo original y el binario deben seguir teniendo los mismos hashes; un cambio exige otra ejecución revisada. Si el objetivo ya terminó, `--resume` devuelve 0 sin evolucionar. Un fallo requiere inspección; el lanzador no convierte un fallo en continuación limpia. Un presupuesto que sólo permita inicializar y guardar el paso cero se considera falta de avance: aumentar el presupuesto y preparar una nueva ejecución.

**En una asignación Slurm usar como máximo un segmento por invocación**, hasta implementar y verificar la integración con el scheduler. Encadenar procesos locales no crea nuevas asignaciones. El límite comunicado debe descontar prólogo, preparación y validaciones del tiempo realmente disponible; el margen debe cubrir checkpoint, diagnósticos, validación y cierre. Los valores75/5/25s son controles locales, no márgenes aceptados para producción de 24h. No hay requeue, señales ni consulta automática al tiempo restante Slurm.

## Archivos para revisar el progreso

- `segments.json`: estado de la cadena, hashes y resultados de cada segmento.
- `segment_0001/run.par`, `run.log`, `check.log`/`verify.log`, `receipt.json`: configuración efectiva, salida y comprobaciones. La primera ejecución no necesita `check.log` de reinicio.
- `data/run_status.dat`: estado del solver y checkpoint publicado más reciente.
- `data/weno5_hlle/`: series y VTK comunes a toda la cadena.

Por ejemplo:

```bash
cat .validation/michel_segmentos_nuevo/data/run_status.dat
python3 -m json.tool .validation/michel_segmentos_nuevo/segments.json
```

No seleccionar temporales `.rst.tmp`. No se elige silenciosamente otra generación si falla la que indica el estado: el error detiene la cadena y conserva la evidencia.

## Estado que conserva V9

El cuerpo V8 se mantiene y se añaden dos registros antes del tiempo/paso/estado físico:

1. Los16contadores enteros de eventos numéricos y las cinco sumas de cambios conservativos por reparación.
2. Bandera de balance radial activo, masa inicial discreta y las dos integrales temporales de flujo radial de Michel.

El cierre ahora es `GRIAL_END_V9`. Lectores V5/V6/V7/V8 permanecen y sus fixtures se verifican; los formatos anteriores a V9 no contienen acumuladores y advierten que comienzan una nueva base diagnóstica. Para continuar una campaña antigua, usar otra carpeta de salidas para los nuevos contadores/balance; no interpretar sus contadores como acumulados históricos recuperados. Los pisos V8/V9 se restauran según el contrato de R01.

La recuperación de primitivas necesaria para arrancar no aumenta los contadores físicos de evolución. El balance radial se restaura antes de integrar nuevos flujos; no se cambia su masa inicial en cada proceso. Si se habilita por primera vez un diagnóstico que antes estaba desactivado, no se reconstruye el historial que no se midió. Mantener los controles diagnósticos fijos para una comparación de continuidad.

## Series y protección de evidencia

Cada serie tiene un último tiempo en caché. Al iniciar se lee su historial una sola vez: se exige orden estrictamente creciente y se rechazan tiempos posteriores al checkpoint. Una muestra inicial de reinicio que ya existe no se vuelve a escribir. Las muestras terminales de una parada sí se conservan: una cadena puede tener más filas que la evolución continua, pero ninguna fila de tiempo repetido y ninguna pérdida del muestreo regular.

GW y acreción ahora serializan tiempo y valores con 17 dígitos significativos; conservan sus columnas. La precisión evita confundir pasos diferentes o comparar un tiempo redondeado con el checkpoint. Una serie antigua de menor precisión puede requerir continuar en otro directorio, conservando el original.

Si las series están adelantadas respecto al checkpoint elegido, se detiene la ejecución; no se recorta ni sobrescribe el historial. Esto puede ocurrir al escoger un checkpoint anterior o después de un corte abrupto que dejó datos posteriores a la última generación. Las figuras VTK existentes tampoco se sobrescriben al reiniciar. La recuperación tras un corte de ese tipo necesita una rama explícita en otra carpeta y revisión del empalme; no está automatizada por `--resume`.

## Evidencia y resultados

Evidencia completa: `.validation/r04_segments_20261024_v1/VERIFIED.json`. Configuraciones reproducibles en `benchmarks/restart/`; scripts, logs, checkpoints, archivos previos, parches y hashes en el dossier local.

| Control | Malla | Tiempo físico final | Pasos | Segmentos | Duración de la cadena | Continua |
|---|---|---:|---:|---:|---:|---:|
| Michel, EF/logarítmico |1024×1×1|40M|32392|4|219.77s|204.07s|
| KHI, cartesiano |200×1×200|2|2667|4|238.13s|215.23s|

Ambas configuraciones usan FP64, WENO5–HLLE, sensor, CFL=0.3 y cuatro hilos. Los tres primeros procesos de cada cadena duraron aproximadamente 70s; el último completó el objetivo. Se fijaron límites 75s, margen 5s y checkpoint periódico 25s. Son pruebas concurrentes locales; las duraciones no son un estudio de escalabilidad.

**Resultado:** campos finales de cinco conservadas, paso/tiempo, contadores y acumuladores de reparación/flujo idénticos por bytes a sus respectivas evoluciones continuas. Todas las muestras comunes comparadas tienen diferencia0, tiempos estrictamente crecientes y los mismos frames VTK (nueve Michel, cinco KHI). Las filas adicionales corresponden a paradas, no a duplicación del empalme. En Michel se conservan M0 y ambas integrales radiales a lo largo de los cuatro procesos.

Controles complementarios:

- 23 regresiones Fortran y dos grupos de integración, optimizados y debug con bounds/traps `invalid,zero,overflow`.
- Michel debug reiniciado desde un checkpoint intermedio: campos, contadores y balance radial idénticos.
- Toro FD 3D 48×24×16 con perturbación en 0.05M, evolucionado hasta 0.2M en dos segmentos: un solo evento, campos idénticos y series PPI/acreción/GW sin tiempos duplicados. Es control de continuidad, no estudio físico prolongado de PPI.
- Lanzador: tres segmentos simulados, cuota y `--resume`, final ya completado, fallo de solver, estado ausente, checkpoint inválido, ausencia de avance, final falso, cambio de ejecutable y protección de directorio existente.
- Escritura realmente interrumpida mediante límite de tamaño de archivo/SIGXFSZ: temporal parcial de 4096 bytes, generación anterior intacta por SHA256 y legible; el parcial se rechaza. Un historial adelantado también se rechaza sin modificar sus bytes.

En KHI se comparan balances y diagnósticos de sus interfaces. Su salida auxiliar `ppi_modes.dat` sólo se usa como control de reproducibilidad; el control físico de modos azimutales corresponde al toro.

**R04 y Control R cerrados localmente. Siguiente tarea: M01**, compilación opcional/API MPI y errores coordinados. No hay todavía comunicación MPI ni prueba de una asignación real de 24h; siguen pendientes calibración de márgenes y V04. El código principal, la referencia y las corridas CADS se mantienen fuera de este incremento.
