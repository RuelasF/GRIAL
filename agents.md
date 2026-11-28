# Contexto de auditoría para GRIAL

GRIAL resuelve hidrodinámica relativista en FP64 sobre una métrica analítica fija. Para auditar MPI y rendimiento, usa el núcleo de `mpi/` y toma el intercambio bloqueante como referencia. El objetivo es que un lector pueda seguir un paso de evolución, entender qué comunica cada rango y justificar cada opción adicional con una medida reproducible.

## Base de Git

Esta guía está en `main`. La auditoría de las mejoras hidro/MPI parte de [`hidro-rendimiento` en el commit publicado `0c0fd13`](https://github.com/RuelasF/GRIAL/tree/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi). Los enlaces a código y manuales apuntan a ese commit. Los SHA de código se verificaron el 27 de noviembre de 2026. Antes de empezar, registra el SHA que vas a auditar y cualquier cambio local.

| Rama | SHA verificado | Papel |
| --- | --- | --- |
| `hidro-rendimiento` | `0c0fd1372aa5463c9adb0121247ba599f11cc10b` | Versión hidro/MPI con las correcciones y optimizaciones recientes |
| `main` | `b95e7acb05133feb33196b189517ff1afaf0c900` | Referencia hidro anterior a esas mejoras, antes de publicar esta guía en `main` |
| `mpi-hidro` | `cef9a8750c239fb115a398f005aaac61350e7501` | Referencia histórica de la publicación MPI |
| `grmhd` | `1d843f93985c825ab22396e7ffd7e508d9c1bb53` | Rama magnética experimental separada |

El ejecutable de la raíz es histórico. Comprueba cualquier hallazgo allí en la implementación de `mpi/` antes de atribuirlo a la ruta distribuida. El trabajo magnético local `grmhd-wip` y los prototipos aislados tienen su propia procedencia. Describe los hallazgos magnéticos aparte y cita la versión concreta.

```bash
git clone --single-branch --branch hidro-rendimiento https://github.com/RuelasF/GRIAL.git
cd GRIAL
git switch --detach 0c0fd1372aa5463c9adb0121247ba599f11cc10b
git rev-parse HEAD
git status --short --branch
```

## Ruta de ejecución

| Qué seguir | Archivos |
| --- | --- |
| Entrada `--distributed`, bucle y eventos de salida | [grial.f90](https://github.com/RuelasF/GRIAL/blob/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi/grial.f90), [parallel_application.f90](https://github.com/RuelasF/GRIAL/blob/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi/parallel_application.f90) |
| Inicio MPI y reparto | [parallel_runtime.f90](https://github.com/RuelasF/GRIAL/blob/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi/parallel_runtime.f90), [domain_decomposition.f90](https://github.com/RuelasF/GRIAL/blob/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi/domain_decomposition.f90), [initialization.f90](https://github.com/RuelasF/GRIAL/blob/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi/initialization.f90) |
| Fantasmas y etapas RK3 | [halo_exchange.f90](https://github.com/RuelasF/GRIAL/blob/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi/halo_exchange.f90), [time_integrator.f90](https://github.com/RuelasF/GRIAL/blob/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi/time_integrator.f90) |
| CFL, RHS y recuperación de primitivas | [parallel_timing.f90](https://github.com/RuelasF/GRIAL/blob/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi/parallel_timing.f90), [evolution.f90](https://github.com/RuelasF/GRIAL/blob/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi/evolution.f90), [equations.f90](https://github.com/RuelasF/GRIAL/blob/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi/equations.f90) |
| Balances, observadores, VTK y reinicio | [parallel_diagnostics.f90](https://github.com/RuelasF/GRIAL/blob/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi/parallel_diagnostics.f90), [parallel_observer.f90](https://github.com/RuelasF/GRIAL/blob/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi/parallel_observer.f90), [parallel_visualization.f90](https://github.com/RuelasF/GRIAL/blob/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi/parallel_visualization.f90), [parallel_checkpoint.f90](https://github.com/RuelasF/GRIAL/blob/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi/parallel_checkpoint.f90) |

En toros 3D, `--grid NR NTHETA` reparte radio y theta. Cada bloque conserva todo el azimut. En 1D y en el jet se reparte el primer eje. El runtime solicita `MPI_THREAD_FUNNELED`. Los hilos OpenMP trabajan en el cálculo local y las llamadas MPI quedan fuera de sus regiones paralelas.

La ruta de referencia usa `MPI_Sendrecv` para intercambiar capas de primitivas antes del primer RHS y después de recuperar primitivas en cada etapa RK. Primero se resuelve la dirección radial, incluidas las fronteras físicas del extremo. El intercambio polar incluye las capas radiales para completar las esquinas. Los rangos comparten el paso CFL y acuerdan los eventos colectivos. El reinicio distribuido valida la generación y las historias asociadas antes de restaurarlas. La guía documenta cómo continuar con otro reparto. [Halos](https://github.com/RuelasF/GRIAL/blob/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi/docs/HALOS_M03.md) · [Entrada distribuida y reinicio](https://github.com/RuelasF/GRIAL/blob/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi/docs/ENTRADA_DISTRIBUIDA_CONTROL_C.md).

## Criterio de simplicidad

- Conserva como base la comunicación bloqueante, el contrato de halos y el núcleo numérico compartido entre OpenMP y MPI. Reutiliza las rutinas que ya existen. Señala duplicación, estados sin lectores, opciones sin uso y capas de abstracción que se puedan retirar con una regresión acotada.
- Justifica una ruta alternativa por un cuello de botella medido. Evalúa primero si se puede resolver con un reparto, afinidad, cadencia de salida o cambio local más simple. Para justificar más estados, buffers o sincronizaciones, muestra qué problema resuelven, en qué condiciones aportan y cuánto cuestan en memoria, mantenimiento y reinicios.
- Mantén las comprobaciones necesarias para evitar divergencias entre rangos, corrupción de checkpoints o aceptación de estados inválidos. Al retirar una colectiva, demuestra qué invariante la vuelve redundante, cómo se conserva la secuencia común y cómo se propaga un fallo de un solo rango.
- Mide el bucle completo a igual trabajo, con la misma malla, métodos, pasos, reparto, afinidad, compilador y cadencia de E/S en cada par. Separa el arranque y la escritura final cuando corresponda. Usa réplicas alternadas, informa mediana y dispersión y registra datos de memoria o comunicación si son el motivo del cambio. El tiempo dentro de MPI incluye espera y las regiones anidadas se solapan.
- Comprueba U/P completos, balances, contadores y reinicio, además del tiempo. Para un cambio que propone conservar la aritmética, busca igualdad bit a bit en el mismo reparto. Al cambiar el reparto, compara campos globales y contadores físicos globales. Si cambia la discretización, el criterio de aceptación debe incluir un error o una convergencia física justificada.

## Evidencia que informa la auditoría

- `0c0fd13` reutiliza invariantes de momento en la recuperación de primitivas y barre el CFL causal sobre el plano radial/polar de la geometría axisimétrica. `7310991` redujo a ese mismo plano las cachés de geometría. Para una malla de 256×128×256, su tamaño calculado pasó de 12.33 GiB a 49.29 MiB. Es una estimación de almacenamiento de la caché. El pico RSS se mide por separado.
- En el piloto hidro CADS hubo 145 ejecuciones. La versión actual redujo el tiempo en 35 de 36 pares frente a la referencia. Las medianas por configuración bajaron 2.8–12.0% en el toro FD y 7.0–23.7% en el jet. La comparación mide el efecto conjunto de las mejoras en segmentos cortos, con una asignación por caso y cantidad de nodos. La evidencia incluye comparaciones de campos y reinicios en Lustre con otro reparto.
- `OVERLAP=1` sigue siendo experimental y su efecto varía con el caso y el reparto. La caché optativa `CFL_CACHE=1` mejoró 14 de los 36 pares frente a la versión actual sin esa opción. Ambos tienen valor predeterminado 0. La agrupación de guardas colectivas redujo llamadas pero aumentó la mediana del bucle de 16.683 a 17.811 s en el A/B local del toro. Sigue vigente la secuencia original de guardas.
- La calibración PPI con el toro FD 96×64×128 y `polar_mesh_h=0.3` midió 15.409 s por 96 pasos en el mejor reparto de un nodo (12 rangos × 3 hilos, `6×2`) y 20.243 s en el mejor de dos nodos. La malla sigue siendo candidata. Crecimiento y saturación de la PPI, convergencia de los observables y rendimiento a 20 000 M quedan por validar. El diagnóstico GW es un proxy cuadrupolar comparativo sobre Kerr fijo.

El clon contiene el código, las regresiones y la [documentación técnica](https://github.com/RuelasF/GRIAL/blob/0c0fd1372aa5463c9adb0121247ba599f11cc10b/mpi/README.md). En la máquina de trabajo, consulta también `docs/CONTEXTO_ACTUAL_GRHD.md` del repositorio original por ruta explícita y los informes `D:/Resultados_GRIAL/RENDIMIENTO_HIDRO_CADS_20261123/INFORME.md` y `D:/Resultados_GRIAL/RENDIMIENTO_PPI_CADS_20261123/INFORME.md` si están disponibles. Los dossiers locales y los resultados de clúster quedan fuera del clon público. Indica cuándo un resultado se leyó de un informe y cuándo se reprodujo durante la auditoría.

## Verificación y entrega de hallazgos

Con GNU Fortran, Make, Python 3 y Open MPI, desde la raíz de un clon de auditoría:

```bash
OMP_NUM_THREADS=1 make -C mpi -j2 ARCH_FLAGS='' test-application
OMP_NUM_THREADS=1 make -C mpi -j2 USE_MPI=1 OVERLAP=0 CFL_CACHE=0 ARCH_FLAGS='' test-application
```

Para cambios numéricos o de comunicación, ejecuta las pruebas pertinentes y las suites `make test` en CPU y MPI. Añade comprobaciones de límites, trampas de punto flotante y una prueba multirrango que atraviese la ruta modificada. `make USE_MPI=1 test-application-mpi` desde `mpi/` amplía la cobertura de aplicación. Limita las pruebas locales a cuatro núcleos físicos. Usa salidas nuevas, conserva el trabajo existente y las fuentes congeladas y evita la opción corta `-u`, incluso agrupada. Una auditoría de solo lectura permite reproducir en copias aisladas. Los cambios de solver, commits, publicaciones y trabajos de clúster requieren una solicitud del usuario.

Relaciona cada hallazgo con la rama, SHA, archivo y función. Describe el mecanismo, las condiciones en que se activa, su impacto y cómo reproducirlo. Para una mejora, incluye la alternativa más simple, la prueba de equivalencia y la medición que justificaría adoptarla. Distingue los defectos reproducidos de las hipótesis y del código histórico.
