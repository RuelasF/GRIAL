# Reinicios y parada por tiempo real — R01–R04

Implementado el 23/oct/2026 exclusivamente en la copia `mpi/`. Sigue siendo un solver OpenMP de un proceso. La ruta principal y las campañas CADS no reciben estos cambios.

## Uso

Añadir al archivo `.par`, además de los parámetros físicos habituales:

```ini
walltime_limit_seconds = 86400
walltime_margin_seconds = 900
walltime_checkpoint_seconds = 1800
```

Las unidades son segundos reales. El ejemplo representa 24 h por proceso, margen de 15 min y checkpoints cada 30 min aproximadamente. Es una configuración de ejemplo: hay que medir el coste del RK, los diagnósticos, la escritura y el cierre en el sistema de destino antes de fijar el margen de producción.

- `walltime_limit_seconds=0` desactiva la parada por tiempo real, por defecto.
- `walltime_margin_seconds=900` por defecto. Con límite activo debe cumplir `0 < margen < límite`.
- `walltime_checkpoint_seconds=0` desactiva la periodicidad real, por defecto. Un valor positivo se atiende después de un RK completo.
- `checkpoint_interval` conserva su significado anterior: intervalo de **tiempo físico**, positivo. Ambos calendarios pueden coexistir; una coincidencia genera un solo checkpoint.

El reloj `SYSTEM_CLOCK` usa enteros de64bits, origen local y gestión del retorno del contador. Cuenta desde el inicio del programa, incluyendo lectura, inicialización, recuperación y E/S. No usa `CPU_TIME`. El límite es **por proceso/segmento**: al reiniciar comienza un reloj nuevo. Los controles se conservan en V8 y pueden sobrescribirse explícitamente al reiniciar.

El presupuesto debe ser el tiempo realmente disponible **al arrancar el programa**. Si una asignación de 24 h ya consumió tiempo en compilación, prólogo o preparación, descontarlo al configurar `walltime_limit_seconds`. Aún no hay consulta automática a Slurm ni lanzador que calcule ese remanente.

## Punto seguro y resultado

La decisión se toma entre pasos completos SSP-RK3. No se modifica `dt` para ajustar un corte de tiempo real; sus tres etapas, condiciones de frontera y recuperación ya terminaron. Antes de comenzar otro paso se comprueba:

```text
elapsed_wall + max_observed_whole_step >= walltime_limit - walltime_margin
```

La reserva adicional usa el paso completo más lento observado, incluida su E/S. La inicialización cuenta contra el presupuesto: si lo agotó, se publica un checkpoint del paso cero. Si se alcanzó el tiempo físico final, prevalece el estado `completed`.

Al parar se atienden los diagnósticos terminales, se publica un checkpoint de ese paso y se ejecuta la limpieza normal. Siempre se guarda también un checkpoint al completar naturalmente el objetivo. Una generación ya escrita en ese paso no se duplica.

`output_folder/run_status.dat` se publica mediante archivo temporal y renombrado:

```text
status=stopped_for_walltime
physical_time=...
step=...
checkpoint=native/checkpoint_weno5_hlle_step_000000123.rst
elapsed_wall_seconds=...
```

Los estados son `running`, `stopped_for_walltime` y `completed`. Cada checkpoint publicado actualiza el estado `running` con su ruta; sólo al cerrar se escribe el resultado terminal. Parada planificada y final físico devuelven código0; se distinguen por `status`. Los fallos de checkpoint/estado devuelven error no nulo; un lanzador futuro deberá comprobar **ambos** código y estado. Slurm puede marcar `COMPLETED` una parada planificada que todavía necesita continuación.

No hay un manejador de señales. No añadir `--signal=USR1@...` a esta versión esperando un checkpoint. Un corte externo no atendido pierde el trabajo posterior a la última generación completa. El margen no garantiza cumplir un plazo ante E/S bloqueada, un primer paso excepcionalmente lento o una inicialización que exceda toda la asignación.

## Formatos V8/V9 y compatibilidad

**Formato actual: V9.** R04 conserva el cuerpo V8 descrito abajo y añade contadores numéricos y acumuladores del balance radial de Michel, con cierre `GRIAL_END_V9`. Las series se continúan sin duplicar el empalme. [Lanzador, contrato y resultados R04](SEGMENTOS_R04.md).

V8 conserva el cuerpo hidrodinámico V7 y añade, antes del tiempo/paso/estado:

1. Pisos efectivos `rho_floor`, `p_floor`, `D_floor`, `tau_floor`.
2. `checkpoint_interval`, `extraction_stride`, `mass_monitor_stride`, `numerical_diagnostics`, `vtk_mapping`.
3. Límite, margen y periodicidad por tiempo real.

Después de las cinco conservadas `[x,y,z,variable]`, incluye un registro `GRIAL_END_V8` con paso, tiempo y dimensiones repetidos. El lector verifica registros, estado finito y concordancia del cierre. `--check` de un reinicio lee también el estado y el cierre, sin crear salidas.

Los pisos primitivos tienen ahora claves `.par` `rho_floor` y `p_floor`, con valores por defecto 1e−10/1e−13. En un arranque frío se derivan `D_floor=rho_floor` y `tau_floor=p_floor*(Gamma/(Gamma-1)-1)`. V8 restaura los cuatro valores efectivos. Cambiar explícitamente `rho_floor` o `p_floor` al reiniciar cambia el piso correspondiente y deriva su piso conservativo; es un cambio deliberado de configuración física/numérica y no una continuación idéntica. Los manifiestos y el resumen registran pisos y controles reales.

Se mantienen lectores anteriores. En V7 y versiones previas no se conocen los pisos históricos: se emite una advertencia y se adoptan 1e−10/1e−13, con los conservativos derivados. Para una campaña antigua con otros pisos se deben indicar explícitamente en el archivo de reinicio; no se infieren de la densidad mínima ni se importa un sidecar de un observador externo. No se ha eliminado la compatibilidad V2–V6, pero las fixtures verificadas en esta etapa son V5/V6/V7 y V8.

## Publicación y recuperación

La escritura usa `checkpoint_<scheme>_step_<N>.rst.tmp` en el mismo directorio. Se comprueba cada escritura y el cierre; sólo después se renombra al `.rst` final. Los checkpoints de pasos anteriores permanecen disponibles. Si falla una escritura o el renombrado, el archivo temporal no se publica y se comunica el error. Los temporales de fallo se conservan para inspección y no deben seleccionarse como reinicio.

Se verifica integridad estructural/completitud; el cierre no es una suma criptográfica de todo el contenido. Tampoco se implementó `fsync` para garantizar supervivencia a una pérdida de energía del almacenamiento. Renombrar en un mismo sistema de archivos evita exponer un archivo final mientras se escribe; no reemplaza una copia de seguridad.

Para continuar manualmente, desde el mismo directorio de ejecución y conservando el directorio de salidas:

```ini
restart = true
restart_file = native/checkpoint_weno5_hlle_step_000000123.rst
walltime_limit_seconds = 86400
walltime_margin_seconds = 900
walltime_checkpoint_seconds = 1800
```

```bash
./grial --check restart.par
OMP_NUM_THREADS=18 ./grial restart.par
```

Usar la ruta real publicada en `run_status.dat`, relativa al directorio desde el que se ejecutó. Si se cambia `output_folder`, preparar su directorio y el subdirectorio del esquema antes de continuar, como exige la ruta de reinicio existente. El lanzador R04 continúa desde la generación validada del estado; no elige silenciosamente otra ante un fallo. No hay reenvío/requeue ni comunicación MPI todavía.

## Verificación y siguiente paso

Evidencia: `.validation/walltime_restart_20261023_v1/VERIFIED.json`, logs `suite_opt_final.log`, `suite_debug_final.log` y `final_reference/REFERENCE_VERIFIED.json`. Fuentes previas, diferencias y hashes permiten revisar/revertir sólo este incremento; no restaurar toda la carpeta ni tocar `reference/`.

-23 regresiones Fortran más un grupo de integración de parada/reinicio, optimizados y debug con límites y traps `invalid,zero,overflow`.
- Ida/vuelta de cuatro pisos efectivos y controles, errores forzados de apertura/publicación sin dañar la generación anterior, rechazo de archivos incompletos y formato/cierre inconsistentes.
- Parada durante evolución y por presupuesto agotado en inicialización; guardado periódico real y validación sin salidas.
- Sod continuo frente a interrumpido/reanudado; cinco controles frente a la referencia congelada y debug; FD/FM con pisos personalizados y reinicio antes/después de la perturbación. Comparaciones de estados, tiempo y paso, no sólo código de salida.

**R04 sigue pendiente:** lanzador de segmentos, continuidad de series sin muestras duplicadas en los empalmes y política de presupuestos/contadores entre procesos. Los contadores numéricos y el balance radial existentes siguen siendo por sesión, como indican sus encabezados; V8 no los convierte en acumulados globales. Las series terminales y el muestreo inicial del reinicio pueden repetir el punto de empalme. Primero cerrar R04/Control R; después iniciar M01. Sin pruebas Slurm de 24 h, MPI o GPU en esta etapa.
