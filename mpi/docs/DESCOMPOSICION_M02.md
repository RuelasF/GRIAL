# M02: descomposición del primer eje de GRIAL

24/oct/2026. Sólo copia hidrodinámica `mpi/`. M02 construye la malla local; el ejecutable aún rechaza evolución y `--check` con varios procesos. No equivale a un solver distribuido ni a una medida de escalabilidad.

## Secuencia aprobada

1. **SRHD1D:** cerrar malla, halos, reconstrucción, CFL y RK con onda advectada y choque, frente a un proceso.
2. **SRHD2D:** repetir equivalencia en los planos x-z y x-y; comenzar con controles cortos de KHI, después jet en coordenadas cilíndricas. Separar métrica plana de coordenadas cartesianas: el jet ya tiene geometría curvilínea.
3. **GRHD:** comprobar primero Michel radial y luego toros FD/FM con reparto radial; fuentes, sondas y diagnósticos deben conservar índices y presupuestos globales.
4. Repartir también theta, resolver polos y diagnósticos distribuidos según T01/D01–D03 antes de producción tridimensional.

La física SRHD ya utiliza la métrica de Minkowski; el diseño de la malla no cambia al pasar a una métrica curva. En M02 se comprueba su caché estática para detectar errores de coordenadas desde ahora. La evolución GRHD distribuida no queda validada con esta prueba.

## Reparto

`domain_decomposition.f90` no depende de MPI ni del estado del fluido. `build_domain(shape,halo,rank,ranks,periodic_x,layout,message)` devuelve una descripción canónica o un mensaje de error; `validate_domain` rechaza descripciones inconsistentes.

Para N celdas globales y P procesos, q=N/P y s=mod(N,P). Los primeros s procesos tienen q+1 celdas; los restantes tienen q. El desplazamiento del proceso p es `p*q+min(p,s)`. El índice global de una celda local i es `offset+i`; las celdas físicas son `1:n_local` y las caras `0:n_local`.

- Sólo se distribuye el primer eje x/r. Los tamaños y/z y sus direcciones degeneradas se conservan.
- `layout%global_shape`, `%local_shape`, `%offset` y `%active_axes` separan dominio global y almacenamiento local.
- `%owns_lower/upper` identifica los extremos globales; `%physical_lower/upper` excluye el eje periódico.
- `%lower_rank/upper_rank` son vecinos lógicos. **-1 significa que no existe vecino; no es una constante MPI ni debe pasarse directamente a una comunicación.** M03 hará esa traducción.
- La periodicidad conecta ambos extremos, incluido vecino propio con un proceso. No aplica todavía condiciones de frontera ni intercambia datos.
- Con varios procesos se exige `floor(N/P)>=max(2,halo)`: cubre el stencil y evita que los kernels actuales interpreten `nx=1` como dimensión inactiva. La ruta serial admite ejes degenerados y mallas pequeñas como antes.
- Se rechazan formas/rangos/profundidades no positivos, rangos inexistentes, bloques insuficientes y extensión de índices que desborde el entero.

## Conexión a la malla existente

`allocate_and_grid()` conserva el arranque serial. La variante `allocate_and_grid(layout)` recibe los parámetros **globales** en `nx,ny,nz` y una profundidad de halo coincidente; valida el layout antes de reservar memoria y después convierte esos tamaños en locales.

`variables.f90` añade `nx_global,ny_global,nz_global,x_global_offset`, direcciones globales activas y propiedad de bordes físicos. Los arreglos mantienen FP64 y `[x,y,z,variable]`; la reserva masiva y la caché usan tamaños locales.

El espaciado se calcula con el tamaño global. `x_min/x_max` y `r_min/r_max` conservan los límites globales. Los centros y caras se calculan con **el mismo índice global y orden de operaciones de la malla serial**. Las coordenadas fantasma se extienden sin envolverlas; M03 comunicará los valores periódicos del estado, no cambiará esa coordenada local.

Las métricas no se modifican: reciben las coordenadas locales correctas. La ubicación de horizonte y sonda se calcula sobre la malla global; `idx_horizon_global/idx_probe_global` y `horizon_owner_rank/probe_owner_rank` identifican cada posición. `idx_horizon/idx_probe` son locales y valen cero fuera del propietario. La inexistencia global de una sonda sigue siendo un error si se pide extracción. Esta propiedad no implementa aún la reducción de acreción distribuida.

No se introduce una opción `.par` de reparto antes de conectar la evolución. Tampoco se modifica el lector de parámetros, el checkpoint o las condiciones físicas para simular un arranque MPI incompleto. Queda pendiente convertir las rutas fatales históricas en fallos coordinados antes de M06.

## Comandos de verificación

Desde `~/Codigo/GRIAL/mpi/`:

```bash
make -j2 ARCH_FLAGS='' all test
make -j2 USE_MPI=1 ARCH_FLAGS='' all test test-mpi test-domain-mpi
```

`test-domain` ejecuta la prueba de malla en un proceso. `test-domain-mpi` lanza la misma prueba con 1, 2 y 4 procesos, usando dos hilos OpenMP por proceso, y prueba rechazos colectivos. Usa MPI local; no envía trabajos Slurm. Un agente con sandbox puede necesitar autorización para los sockets locales.

La prueba reserva primero una malla completa pequeña como referencia y después la malla del bloque que corresponde al proceso. Esa reserva completa es **sólo del test**: la rutina de producción con layout reserva únicamente el bloque local.

## Criterios y evidencia

En `.validation/m02_domain_20261024_v1/` se conservan fuentes anteriores, parches, control previo, manifiestos, logs y resultados. `VERIFIED.json` registra el cierre local:25regresiones Fortran+dos grupos de integración en cada uno de cuatro builds (sinMPI/MPI y optimizado/debug); matrices M01 de23controles y M02 de6controles en ambos buildsMPI;20comparaciones de evolución de un proceso frente al ejecutable previo (Sod/KHI/Michel/FD/FM) exactas también en contadores/presupuestos conunhilo. Los builds debug usan `-O0 -fcheck=all -ffpe-trap=invalid,zero,overflow`. La referencia inmutable y las fuentes/binario del proyecto superior se comprueban porSHA256.

Nueve configuraciones por cada profundidad1/2/3, 27 escenarios por proceso:

- Sod1D, KHI2D x-z y dato plano2D x-y.
- Minkowski cilíndrica2D.
- Eddington–Finkelstein radial lineal, Kerr–Schild radial logarítmica a=0 y caché3D a=0.9 con polos.
- Horizonte excisado y horizonte más allá del dominio.

Se comprueban cobertura única y continuidad de bloques1–5, reparto desigual N=37, mínimo admisible, vecinos/periodicidad, propietarios, dimensiones activas, coordenadas de centros/fantasmas/caras y caché completa en centros y tres familias de caras. Las coordenadas y la caché exigen igualdad exacta con la porción serial. Se verifica también el dato inicial y U(P), sin avanzar el tiempo.

**Redondeo previo de KHI:** al variar la longitud de la tira, las rutinas de inicialización sin modificar alternan `sin/exp` vectorizadas y escalares. Se reprodujo con fuentes anteriores y reserva manual, sin el módulo de dominio (`baseline_initial/pre_m02.log`; el binario enlaza ambas variantes libm). La discrepancia observada es8.673617379884035e-19 en velocidad transversal y1.3877787807814457e-17 en conservadas:1 y2ULP. La nueva prueba exige igualdad exacta en densidad/presión/velocidades longitudinales y acota el redondeo transversal y propagado en U(P) a4/16ULP por elemento; registra máximos por rango. Los otros datos iniciales y toda la geometría exigen igualdad exacta. No se cambian algoritmos físicos ni tolerancias de las regresiones anteriores. Las ejecuciones nuevas de un proceso se comparan además exactamente con el ejecutable pre-M02 de las mismas flags.

Las sumas atómicas OpenMP y las diferencias O0/O3 documentadas en M01 se mantienen separadas: las comparaciones de evolución usan un hilo y referencias con flags equivalentes. No se interpreta este control como identidad de evoluciones turbulentas largas.

## Reversión y siguiente etapa

`before/` conserva los archivos modificados y `changes.patch` el incremento. Para revertir, revisar cambios posteriores y restaurar sólo este incremento en `mpi/`; no tocar `reference/`, el código superior ni fuentes congeladas de CADS.

**Siguiente: M03**, intercambio de halos y condiciones físicas en los extremos globales. Después M04 interfaces/sensor, M05 CFL y M06 RK. No hay aún resultados de evolución SRHD/GRHD distribuida ni benchmark multinodo; la campaña de reinicios CADS mantiene su aceptación independiente.
