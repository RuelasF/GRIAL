# T01 — dominio r×theta con phi completo

24/oct/2026. Implementación verificada localmente en CPU, dentro de la copia hidrodinámica `mpi/`. Cierre: `.validation/t01_rtheta_20261024_v1/VERIFIED.json`. La entrada multiproceso del programa principal conserva su guarda hasta completar diagnósticos y E/S distribuidos.

## Contrato del reparto

`build_domain` admite los argumentos opcionales `process_shape=[Pr,Ptheta]` y `periodic_y`. Las llamadas anteriores conservan `[ranks,1]`. El orden es `rank=cr+Pr*ctheta`: radio varía más rápido. El tercer eje permanece completo en cada bloque.

Cada eje reparte el cociente entero y asigna una celda adicional a los primeros bloques del resto. El bloque mínimo de una dirección repartida debe contener al menos `max(2,nghost)` celdas; una dirección activa nunca se convierte en una dimensión local de una celda. Esto conserva la interpretación de los selectores dimensionales existentes. El tamaño global, las direcciones activas, las coordenadas de las caras y la separación de malla son comunes.

`x_global_offset` y `y_global_offset` convierten índices locales en globales. `process_coords`, `lower_neighbor(2)` y `upper_neighbor(2)` describen ambos ejes; `lower_rank/upper_rank` conservan los vecinos radiales de la API anterior. Una interfaz MPI no activa las protecciones del sensor propias del borde global: las rutas CPU por tira y por cara reciben también el índice polar global.

`domain_owner(layout,i,j)` devuelve el propietario único de una celda global. Sin `j`, devuelve el propietario radial dentro de la columna polar del bloque consultante. Por ello el horizonte y la sonda radial tienen un propietario por columna polar. La integración de esas contribuciones pertenece a D02; estos índices no constituyen un diagnóstico global de acreción.

## Fronteras e intercambio

`fill_primitive_halos` ejecuta, cuando se reparte theta:

1. Intercambio radial sobre las celdas físicas transversales.
2. Fronteras radiales físicas en sus propietarios, con las fronteras polares suprimidas.
3. Intercambio polar que incluye los halos radiales ya válidos.
4. Fronteras polares físicas en sus propietarios y periodicidad azimutal.

Las esquinas viajan en las franjas polares; no hay mensajes diagonales. Las máscaras opcionales `fill_lower_y/fill_upper_y` de `set_boundary_conditions` afectan al segundo eje. El eje phi conserva su llenado local. En FD/FM 3D se mantiene el cruce `k -> k+Nz/2` y las paridades existentes de las velocidades, sin introducir otra convención. Se exige `Nz` par. Los cortes sagitales mantienen su reflexión polar.

La rutina de bajo nivel `exchange_primitive_halos` conserva el intercambio radial de M03; se debe usar `fill_primitive_halos` para completar el llenado en ambos ejes y el tratamiento de las esquinas.

Las franjas persistentes son FP64: dos buffers radiales de `g*ny*nz*neq` y, al repartir theta, dos de `(nx+2g)*g*nz*neq`. Cada eje usa dos `MPI_Sendrecv`; se comprueban rango de los contadores, configuración común y tamaño recibido. Los bloques desiguales conservan tamaños de mensaje compatibles entre vecinos de cada eje. No se añaden reducciones globales por intercambio normal. El RK de M06 sigue llenando halos antes del primer RHS y después de cada recuperación.

Se conservan las tres capas utilizadas por MP5/WENO5/sensor. Las comparaciones de transporte cubren las capas efectivamente rellenadas `1-g:N+g`; el plano adicional de reserva `-g`, ajeno al stencil físico, no se trata como dato intercambiado.

## Casos y límites

| Reparto polar | Condiciones |
| --- | --- |
| `Conv`/`6`, `KHI`/`4` | Periodicidad en los dos primeros ejes |
| `FontDaigneSag`/`11`, `FishMoncSag`/`14` | Radio físico y reflexión polar |
| `FontDaigne3D`/`12`, `FishMonc3D`/`15` | Radio físico, cruce polar y phi completo con `Nz` par |

Otros casos conservan el reparto radial anterior y se rechazan si se solicita reparto polar. No hay selección automática del grid de procesos ni nuevos parámetros de producción. El transporte es CPU; la ruta por cara se comprueba en el host, sin ejecutar GPU.

## Reproducir los controles

Desde `~/Codigo/GRIAL/mpi/`:

```bash
make -j2 ARCH_FLAGS='' test-rtheta
make -j2 USE_MPI=1 ARCH_FLAGS='' test-rtheta-mpi
```

El ejecutable de prueba admite `Pr Ptheta report`, por ejemplo:

```bash
OMP_NUM_THREADS=2 mpirun --oversubscribe --bind-to none -n 8 \
  ./.build/mpi/tests/test_rtheta 2 4 report
```

Ejecutar el ejemplo en un directorio vacío, indicando una ruta absoluta al ejecutable si se cambia de directorio: los informes por rango usan `status='new'` para proteger resultados anteriores. El harness Python crea directorios separados, impone tiempos máximos y conserva los rechazos esperados.

El driver compara contra un dominio completo de 37×19×8, o 37×19×1 para los cortes sagitales. Incluye onda periódica tridimensional, KHI, FD/FM sagitales y FD/FM 3D con a=0/0.9; cinco reconstructores con sensor y HLLE/HLLC. Comprueba geometría en centros/caras, franjas codificadas con dos épocas y capas 1/2/3, esquinas, cruce polar, RHS por tira/cara, primitivas/conservadas evolucionadas, secuencia completa de dt, positividad y cuentas locales de recuperación. Para aislar RK/halos, las condiciones físicas iniciales de la evolución se copian por índices desde la referencia completa, como en M06. No se cierra D03 con este control.

Los tiempos físicos son 0.03 en los casos periódicos y 0.1M en los toros: son pilotos de equivalencia del reparto. La conservación de las cinco variables se comprueba en la onda periódica con umbral 1e-11. No se exige conservación sin flujos a un toro abierto, ni se usan estos pilotos como validación de equilibrio prolongado o PPI.

Cada rango vuelve a calcular el dominio completo de referencia antes de evolucionar su bloque. El coste total de esta comprobación crece con el número de procesos y no representa el rendimiento del solver distribuido. Ejecutar las matrices locales una por vez evita que distintas pruebas compitan por los mismos núcleos.

## Resultados verificados

En cada build, grids 1×1,2×1,1×2,4×1,1×4,2×2,4×2,2×4; dos hilos OpenMP por rango. 80 trayectorias por grid y 510 pasos completos por rango sumando los 80 casos; 1280 comparaciones opt/debug. Campos primitivos/conservados, dt, coordenadas, métricas y halos coinciden exactamente con su referencia completa. RHS por tira y por cara en el host: error escalado cero, con umbral 2e-13. Onda periódica sin reparaciones y máximo error relativo de conservación de las cinco variables 6.29e-15, frente al umbral 1e-11. Los recuentos físicos de recuperaciones y reparaciones coinciden exactamente al agregarlos fuera del solver para los ocho grids de cada build.

Dieciséis rechazos por build: bloque polar pequeño, grid incompatible, periodicidad incorrecta, configuración entre rangos, caso no soportado, Nz impar, layout alterado y cambio del offset polar tras inicializar halos. Sin bloqueos. Las suites completas pasan 30 programas Fortran y dos grupos de integración en cuatro builds serial/MPI opt/debug; matrices previas runtime 23/dominio 6/halos 23/interfaces 7/tiempo 47/RK 19 en ambos modos. Veinte evoluciones del programa habitual conservan campos, contadores y presupuesto radial frente a pre-T01.

Se conservan los intentos iniciales. Un fixture inválido cambió Ny global sin copiarlo a los parámetros de asignación y se detenía antes del rechazo colectivo previsto; se corrige sólo esa rama. Los controles positivos no cambian y su identidad de fuente se comprueba antes de reutilizar sus resultados. Se repiten los dieciséis rechazos en ambos builds. Una matriz RK anterior agotó su temporizador bajo concurrencia de launchers; su repetición con afinidad solicitada libre y OpenMP pasivo pasa con el mismo límite, sin relajar tolerancias numéricas.

## Evidencia y reversión

El directorio `.validation/t01_rtheta_20261024_v1/` conserva las fuentes anteriores, un ejecutable serial previo, fuentes debug, logs y resultados separados, parches antes/después y manifiestos SHA256. La referencia hidrodinámica inmutable y los archivos protegidos del proyecto principal se verifican por hash. Revertir requiere revisar el parche contra el estado actual; no sobrescribir avances posteriores ni campañas congeladas.

El siguiente paso previsto es D01: balances con flujos ponderados por RK, normas, extremos y contadores globales. Siguen pendientes D02/D03, checkpoints/reinicios/E/S distribuidos y la validación de producción en clúster. No se mide escalabilidad con estas pruebas locales.

Las recuperaciones y reparaciones de celdas se pueden sumar entre bloques sin contar caras compartidas. Las llamadas a reconstrucción/HLLC incluyen evaluaciones adicionales en interfaces MPI y deben conservarse como medidas de trabajo, no interpretarse como número de celdas físicas. La auditoría de estos informes es externa al solver: no implementa todavía D01. Los toros gruesos contienen reparaciones de atmósfera; la coincidencia del reparto no supone que esas reparaciones desaparezcan ni certifica su equilibrio prolongado.
