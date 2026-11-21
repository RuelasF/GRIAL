# Comparación MPI de choques SRHD 1D hasta t=0.4

24/oct/2026. M06B, antes de ampliar la descomposición al ángulo polar. La comparación MPI está verificada: 128 trayectorias principales, 32 referencias seriales y 96 comparaciones con MPI 1/2/4, sin diferencias en campos ni trazas temporales. Se utilizan los kernels y el RK de M06, sin cambios numéricos.

## Figuras para revisar

Los perfiles principales incluyen la solución exacta, las tres cantidades de procesos y una fila inferior con las diferencias respecto a MPI 1. Las curvas coincidentes se distinguen con línea continua, discontinua y marcadores huecos; no se desplazan artificialmente.

- [Sod: WENO5–HLLE y diferencias MPI](../.validation/m06b_shocks_t04_20261024_v1/comparison_v2/sod_weno5_hlle_mpi.png).
- [Choque fuerte: WENO5–HLLE y diferencias MPI](../.validation/m06b_shocks_t04_20261024_v1/comparison_v2/strong_weno5_hlle_mpi.png).
- [Sod, cinco reconstructores con HLLE](../.validation/m06b_shocks_t04_20261024_v1/comparison_v2/sod_hlle_reconstructores.png) y [con HLLC](../.validation/m06b_shocks_t04_20261024_v1/comparison_v2/sod_hllc_reconstructores.png).
- [Choque fuerte, cinco reconstructores con HLLE](../.validation/m06b_shocks_t04_20261024_v1/comparison_v2/strong_hlle_reconstructores.png) y [con HLLC](../.validation/m06b_shocks_t04_20261024_v1/comparison_v2/strong_hllc_reconstructores.png).

Las dos figuras principales también están en PDF. En los paneles inferiores, una línea en cero acompañada de «máx |diferencia| = 0» significa igualdad exacta en los datos FP64 reconstruidos, no que las soluciones se hayan sustituido por una curva común.

## Configuración física y numérica

Minkowski cartesiano, dominio [0,1], interfaz inicial x=0.5, velocidad inicial nula y fronteras de salida. Son los tiempos y parámetros canónicos de la validación SRHD previa. Sod utiliza Gamma=1.4; el fixture corto de M06 había usado Gamma=5/3 y no se identifica como la misma corrida.

| Caso | Gamma | rho izquierda/derecha | presión izquierda/derecha | Tiempo final |
|---|---:|---:|---:|---:|
| Sod | 1.4 | 1 / 0.125 | 1 / 0.1 | 0.4 |
| Choque fuerte | 5/3 | 1 / 1 | 1000 / 0.01 | 0.4 |

- N=400 para Godunov, TVD–MC, WENO3, MP5 y WENO5, con HLLE/HLLC y sensor activo.
- N=401 para MP5/WENO5–HLLE, con el mismo dominio y tiempo. Se prueba reparto desigual: 201/200 celdas con dos procesos y 101/100/100/100 con cuatro.
- CFL=0.3; halo de tres capas; pisos rho=1e−10, p=1e−13, D=1e−10 y tau=3e−13; umbrales del sensor por defecto.
- Un hilo OpenMP por proceso para esta comparación 1D. Las cantidades 1/2/4 de las figuras son procesos MPI, no hilos de un proceso.
- Debug con bounds y traps invalid/zero/overflow: MP5/WENO5–HLLE, ambos casos y ambas mallas, también con referencia serial y MPI 1/2/4.

En N=400 la discontinuidad coincide con una cara. En N=401 atraviesa la celda central; el driver integra su fracción izquierda/derecha para inicializar los promedios conservados. Como la velocidad inicial es cero y la EoS Gamma es común, el promedio lineal de rho y p produce exactamente esos promedios de D y tau. El inicializador de producción no se modifica.

Cada ejecución completa 534 pasos con N=400 y 535 con N=401. El último dt se recorta a t=0.4; no se recortan pasos para generar salidas intermedias. Se usa la misma rutina `advance_rk_step` que el programa principal y se coordina la época con M05.

Las ondas exactas de Sod están en x≈0.27689, 0.44366, 0.67041 y 0.78956; las del choque fuerte en x≈0.17347, 0.76725, 0.88416 y 0.89472. Las estructuras están desarrolladas y todavía dentro del dominio físico. Extender arbitrariamente este problema hasta t=1 o 2 haría que las ondas alcanzaran los bordes y exigiría cambiar el dominio o la referencia.

## Resultados MPI y referencia exacta

La comparación reúne las celdas por índice global, comprueba cobertura sin huecos/duplicados y exige igualdad por bytes de coordenadas, cinco primitivas y cinco conservadas respecto a una ejecución serial del mismo build. También exige la misma secuencia completa de dt/tiempos en todos los rangos y en la referencia. Los archivos decimales conservan 18 cifras significativas para recuperar FP64.

Resultado: error máximo MPI–serial **0**, tiempos y campos idénticos, estados finitos/positivos/subluminales y **cero reparaciones de primitivas/conservadas** en las 96 comparaciones. El contador de recuperación es exactamente 3×N×pasos; sólo cuenta celdas físicas. Los contadores de reconstrucción/Riemann sumados representan llamadas, incluidas las caras internas evaluadas por ambos vecinos, y no se presentan como caras físicas únicas.

Como control independiente del driver, el programa habitual `grial` evoluciona cuatro casos completos N=400: Sod/choque fuerte con MP5/WENO5–HLLE. Sus conservadas finales en checkpoint coinciden exactamente con las del driver. [Recibo](../.validation/m06b_shocks_t04_20261024_v1/MAIN_PROFILE_VERIFIED.json).

Se reutiliza la [referencia exacta independiente](../../docs/srhd_1d/exact_riemann.py). Se integran las conservadas exactas dentro de cada celda con cuadratura de orden 12, separando en los frentes. Las normas se calculan directamente contra `up`, sin reconvertir el perfil primitivo. `REFERENCE.json` registra parámetros, posiciones y SHA256 del generador. Los perfiles primitivos de las figuras muestran la solución exacta puntual.

Ejemplo de L1 frente a medias conservadas exactas, WENO5–HLLE, N=400; estos errores son comunes a MPI 1/2/4:

| Caso | L1(D) | L1(tau) | L1(Sx) |
|---|---:|---:|---:|
| Sod | 0.00167309 | 0.00299929 | 0.00240212 |
| Choque fuerte | 0.465762 | 6.466985 | 5.949970 |

La igualdad MPI no implica igualdad con la solución exacta. En el choque fuerte, el intervalo entre contacto y choque mide aproximadamente 0.01056: sólo 4.2 celdas a N=400. El pico de densidad queda ensanchado y subestimado en todas las particiones; WENO5–HLLE da rho_max≈6.78654. Las diferencias de discretización y de esquema aparecen en las figuras y en `errors_exact.tsv`. La malla de 401 celdas controla el reparto desigual; no constituye por sí sola un estudio de convergencia.

## Diagnóstico del balance con extremos constantes

Mientras los extremos numéricos permanezcan en sus estados iniciales, la integral ideal satisface

\[
\langle D\rangle=(\rho_L+\rho_R)/2,\quad
\langle\tau\rangle=(p_L+p_R)/[2(\Gamma-1)],\quad
\langle S_x\rangle=t(p_L-p_R),\quad
\langle S_y\rangle=\langle S_z\rangle=0.
\]

Se conserva el criterio escalado 1e−11 para este diagnóstico. Seis comparaciones lo exceden: las tres particiones del choque fuerte Godunov–HLLE y Godunov–HLLC a N=400. El máximo es 5.54368e−11, idéntico entre particiones. Estas filas quedan marcadas en el informe; no se concede aprobación del balance de extremos constantes a esos controles.

En Godunov–HLLE, el extremo izquierdo final tiene p≈999.999995934 y vx≈1.99e−9. Su flujo energético ya es ≈4.98e−6, frente al cero inicial. La cola numérica de primer orden altera la hipótesis de extremos constantes antes de que llegue el frente analítico. El primer análisis se detuvo por este criterio; se conserva en `comparison/` y `comparison.log`. El análisis final sigue comprobando el mismo umbral y registra por separado el cambio del flujo físico en los extremos.

Control adicional con el mismo Godunov–HLLE serial hasta t=0.4:

| N | Residuo escalado respecto a extremos constantes | Estado izquierdo p, vx | Reparaciones |
|---:|---:|---|---:|
| 400 | 5.54368e−11 | 999.999995934, 1.99e−9 | 0 |
| 800 | 2.84220e−15 | 1000, 0 | 0 |

[Evidencia del refinamiento](../.validation/m06b_shocks_t04_20261024_v1/BOUNDARY_TAIL_REFINEMENT.json). Esta comprobación apoya la atribución del residuo a la cola numérica en el borde. No sustituye integrar los flujos numéricos de cada etapa RK; ese presupuesto general pertenece a D01. La equivalencia MPI se verifica independientemente y da error cero, incluidos los controles Godunov señalados.

## Reproducción, evidencia y alcance

`analysis/mpi_shock_profiles.f90` es un driver de análisis con presets para estos dos choques. Produce perfiles, trazas y metadatos por rango en carpetas nuevas; no utiliza los escritores de producción ni admite reinicios. `tools/compare_shocks_mpi.py` ejecuta la matriz y genera tablas y figuras. No se ha retirado la guarda multiproceso de `grial`.

```bash
cd ~/Codigo/GRIAL/mpi
make -j2 ARCH_FLAGS='' analysis/mpi_shock_profiles
make -j2 USE_MPI=1 ARCH_FLAGS='' analysis/mpi_shock_profiles
python3 tools/compare_shocks_mpi.py --help
```

El comando completo y las compilaciones debug se conservan con sus logs en `.validation/m06b_shocks_t04_20261024_v1/`. El output debe ser una carpeta nueva. Se necesitan numpy, scipy y matplotlib, disponibles en el entorno científico utilizado. El driver MPI está en `.build/mpi/analysis/mpi_shock_profiles`; no se genera un alias de ese archivo en `analysis/`.

[Manifiesto del cierre](../.validation/m06b_shocks_t04_20261024_v1/VERIFIED.json), [comparaciones](../.validation/m06b_shocks_t04_20261024_v1/comparison_v2/VERIFIED.json), [tabla MPI](../.validation/m06b_shocks_t04_20261024_v1/comparison_v2/comparisons.tsv) y [normas exactas](../.validation/m06b_shocks_t04_20261024_v1/comparison_v2/errors_exact.tsv).

Se repitió `make test` sin MPI, optimizado y debug: 29 programas Fortran y dos grupos de integración por build. Los drivers MPI de esta campaña se ejecutaron en optimizado/debug con sockets locales fuera del sandbox. Fuentes Fortran previas, tests existentes, referencia inmutable y archivos protegidos del padre permanecen intactos. Se conservan copias previas, parches, manifiestos y los intentos iniciales; el primer control del programa habitual usó el alias inválido `strong` y se corrigió a `strong_shock` en un directorio nuevo.

Este cierre es de equivalencia y comparación de choques SRHD 1D, con el diagnóstico de borde señalado. No mide escalabilidad ni valida GPU, multinodo, observadores globales, E/S o reinicios MPI de producción. No se consultaron ni modificaron CADS, trabajos remotos, Git, D: o el desarrollo magnético. T01 continúa pendiente; se entregan las figuras antes de comenzar el reparto radial/polar.
