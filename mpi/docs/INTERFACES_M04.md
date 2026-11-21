# M04: reconstrucción y flujos en interfaces MPI

24/oct/2026. Sólo copia hidrodinámica `mpi/`, CPU. M04 usa la posición global para decidir las protecciones geométricas de la reconstrucción y verifica estados, flujos y RHS estáticos en bloques. La evolución multiproceso del ejecutable sigue bloqueada hasta M06.

## Problema y corrección

Los kernels identificaban la cara radial interna por `i<=0` y la franja de protección de MP5 por `i<=3`. Al partir la malla, cada bloque tiene su propia cara local 0: en geometría curva se activaban protecciones sobre interfaces internas. Las primitivas podían estar bien comunicadas y aun así reconstruirse con métodos distintos a ambos lados de la misma cara.

Ahora los datos del stencil conservan sus índices locales, mientras las decisiones geométricas usan `global_face=face_offset+local_face`. Por ejemplo, para N=37 y cuatro bloques, los offsets son 0/10/19/28: la cara local 0 del segundo bloque es la cara global 10 y no es el extremo interno del dominio.

No se cambian los polinomios, limitadores, umbrales, fórmula del sensor, unión de sus seis muestras, recuperación de primitivas, fuentes ni solvers de Riemann. La ruta CPU serial es la referencia.

## Interfaces de código

La reconstrucción por tira admite dos argumentos opcionales:

```fortran
call reconstruct_1d_core(line,n_local,q_left,q_right,direction,face_offset,global_nodes)
```

Sin ellos conserva origen 0 y tamaño global igual al tamaño de la tira, como en la ruta serial previa. Con ellos comprueba que el bloque esté dentro del dominio global. Las posiciones locales siguen determinando los accesos a memoria; sólo cambian las condiciones de selección por geometría. La auditoría de métodos recibe el mismo origen y tamaño, para no registrar un falso fallback en el inicio de cada bloque.

`reconstruct_face_state` mantiene su firma. Sus argumentos `face_index,nodes` se interpretan como posición y tamaño **globales**, puesto que su stencil ya contiene todos los valores necesarios y no se indexa con esa posición.

En `evolution.f90`:

- El barrido CPU radial pasa `x_global_offset` y `max(nx,nx_global)`. El máximo conserva el uso de fixtures seriales que reservan arreglos directamente y dejan `nx_global=0`; una malla real de M02 ya tiene metadatos consistentes, comprobados en M03.
- La ruta por cara captura origen y tamaño radial en escalares anfitriones y los pasa mediante `firstprivate` al bucle OpenACC. Sus llamadas usan la posición global.
- Theta y phi siguen completos en cada proceso, por lo que sus metadatos actuales son0 y el tamaño completo del eje. Las pruebas de ventanas polares verifican la API genérica; no implementan distribución de theta.

## Concordancia entre tira y cara

La revisión reprodujo dos diferencias previas en geometría curva con sensor apagado:

1. TVD por cara imponía Godunov en la cara radial 0, mientras la ruta CPU sin sensor conservaba TVD.
2. MP5 por cara imponía WENO3 en las franjas geométricas aun sin sensor; CPU conservaba MP5 allí, salvo su protección Godunov en la cara radial 0.

Se alineó la ruta por cara con CPU. Con sensor activo se mantienen las protecciones existentes: Godunov en la cara radial interna; para MP5, WENO3 en las primeras tres caras radiales y en las franjas polares especificadas por la política serial. Con sensor apagado se conserva exactamente la política CPU anterior.

El control previo guardado en `PREEXISTING_MISMATCH.log` encuentra diferencias de hasta 1.6832351256068634e-4 por tratar el inicio de un bloque como frontera, y 1.52587890625e-5 entre cara/TVD y tira en el ejemplo curvo. El mismo control corregido, con los metadatos globales y las políticas alineadas, da 0 en `POST_POLICY.log`.

La corrección por cara se realiza sólo en esta copia MPI. Se prueba su ejecución anfitriona; no se ha compilado ni ejecutado en GPU durante M04 y no hay GPU distribuida.

## Comandos

Desde `~/Codigo/GRIAL/mpi/`:

```bash
make -j2 ARCH_FLAGS='' all test
make -j2 USE_MPI=1 ARCH_FLAGS='' all test test-mpi test-domain-mpi test-halos-mpi test-interfaces-mpi
```

`test-interfaces` ejecuta el control de un proceso. `test-interfaces-mpi` lanza procesos reales locales 1/2/4, dos hilos OpenMP por proceso y rechazos de metadatos incorrectos. No envía trabajos Slurm. Un agente con sandbox necesita ejecución fuera de éste cuando los sockets MPI locales están bloqueados.

Los errores de argumentos de la rutina de reconstrucción usan `error stop`, como otros kernels heredados. El contexto de producción se valida colectivamente en M03. Los controles 2/4 con un solo rango inválido terminan sin bloqueo con el lanzador OpenMPI local; la revisión general de rutas fatales/coordinación sigue siendo requisito antes de M06.

## Matriz y criterios

`tests/test_mpi_interfaces.f90` comprueba:

- 1176 ventanas por proceso: dos geometrías de la política (plana/curva), dos juegos de umbrales A/B, 14 configuraciones de método/limitador/sensor, tres ejes y siete offsets, incluidos los dos extremos globales y ventanas interiores. Se comprueban estados y contadores de selección del método.
- 25 fixtures de malla: cinco geometrías por cinco perfiles. Geometrías: SRHD 1D cartesiana, KHI2D periódico, jet2D cilíndrico, corte meridional FD en KS con a=0 y malla3D FM en KS con a=0.9. Las etiquetas FD/FM seleccionan sus fronteras; los campos de estos fixtures son sintéticos, no datos iniciales de equilibrio de toros.
- Perfiles suave, salto débil, medio, fuerte y franja de atmósfera. Los saltos sin sensor prueban equivalencia numérica; no se aceptan como controles de estabilidad, en particular MP5 sin sensor.
- 14 configuraciones: Godunov, TVD con minmod/superbee/MC, WENO3, MP5 y WENO5, con y sin sensor. Cada una usa HLLE y HLLC.
- Reparto radial desigual N=37 y halos de tres capas comunicados por M03. Se comparan todos los estados y flujos radiales con la malla completa y se intercambian directamente los resultados de las caras compartidas para comprobar concordancia entre vecinos.
- El RHS local completo incluye divergencias transversales y fuentes geométricas y se compara con la porción serial. También se compara la ruta anfitriona por cara, con memoria local.

Las comparaciones MPI de estados, flujos, RHS y caras compartidas exigen error 0. La comparación entre implementaciones por tira/cara conserva los criterios de sus regresiones previas, con escala `max(1,norma)`: 5e-15 en estados y 2e-13 en flujos/RHS. Los máximos observados en ambas compilaciones MPI son 0 también en estas comparaciones. No se amplían tolerancias de pruebas existentes.

La cuenta local de métodos se verifica por ventana; las caras MPI aparecen en ambos procesos. La reducción y eliminación de doble conteo global siguen pendientes de D01, al igual que presupuestos y observables distribuidos.

## Evidencia y conservación de la referencia

En `.validation/m04_interfaces_20261024_v1/`:

- `VERIFIED.json`: cierre local y matrices aceptadas; 27 regresiones Fortran+dos grupos de integración en cada uno de cuatro builds sinMPI/MPI optimizado/debug. Debug con límites y traps invalid/zero/overflow.
- Matrices MPI de infraestructura 23, dominio 6, halos 23 e interfaces 7 en cada buildMPI. Se conservan recibos por proceso porque stdout puede intercalar líneas.
- `CPU_DEFAULT_VERIFIED.json`: estados reconstruidos en todas las direcciones activas y 700 RHS completos porbytes frente a una compilación de `before/reconstruction.f90` y `before/evolution.f90`, con iguales flags y un hilo. Incluye 14 configuraciones, cinco geometrías, cinco perfiles y ambos solvers.
- `HYDRO_VERIFIED.json`: 20 comparaciones de evolución de un proceso en Sod/KHI/Michel/FD/FM frente al ejecutable pre-M04 de las mismas flags, exactas en campos/contadores/presupuestos. Se mantienen separadas las diferencias O0/O3 y sumas atómicas OpenMP ya documentadas.
- `before/`, `after/`, `changes.patch`, manifiestos de fuentes/ejecutables y logs de compilación. Referencia inmutable y archivos protegidos del proyecto superior comprobados por SHA256.

No se modifica el ejecutable principal de la carpeta superior ni las campañas congeladas de CADS. Sus estados anteriores mantienen fecha y su aceptación sigue siendo independiente.

## Siguiente y reversión

**M05:** máximos direccionales globales antes de formar el CFL, sincronización de eventos y decisiones temporales. Después M06 conectará halos en cada etapa RK y comparará evoluciones SRHD 1D antes de ampliar a 2D y GRHD. M04 verifica un RHS estático; no equivale a una evolución distribuida ni a un benchmark multinodo.

Para revertir, revisar cambios posteriores y aplicar sólo este incremento desde `before/` y `changes.patch` en `mpi/`. No restaurar globalmente el repositorio ni editar `reference/` o las fuentes congeladas de CADS.
