# D02: observables distribuidos y finalización en rango 0

Estado: verificado localmente en optimizado y debug, con y sin MPI. Evidencia propia: `.validation/d02_observables_20261024_v1/`. El cierre D01 anterior conserva sus fuentes, resultados y manifiestos.

## Reparto del cálculo

Cada proceso integra exclusivamente sus celdas físicas, sin halos. `MPI_Reduce` reúne los pocos coeficientes/componentes de salida en rango 0; allí se calculan amplitudes, fases y polarizaciones y se escriben las series. No se reúne una copia de la malla global. Todos los procesos deben participar en las llamadas colectivas, aunque sólo rango 0 use el resultado final.

La decisión responde a la preferencia del usuario de terminar los diagnósticos en rango 0. No se ha medido una mejora de rendimiento frente a Allreduce, ni se ha confirmado la rutina concreta de OllinSphere. CFL, solicitudes de parada y valores intermedios necesarios en todos los procesos conservan sus comunicaciones globales.

Se reutilizan las integrales seriales de `output.f90`. Las llamadas habituales sin los nuevos argumentos opcionales conservan su comportamiento y formatos. No se añade otro solver ni se cambian los flujos, las velocidades de onda, el recuperador o el RK.

## PPI

En cada bloque se integra por sector azimutal la masa densificada `up(eq_de)*dx*dy*dz` sobre la máscara histórica `rho > 1000*rho_floor`. Se forman las contribuciones a

    C0 = integral D dV_logical,
    Re(Cm) = integral D cos(m phi) dV_logical,
    Im(Cm) = -integral D sin(m phi) dV_logical, m=1..4.

Se reducen **C0, Re(Cm), Im(Cm)**. Sólo después rango 0 obtiene `|Cm|/C0` y la fase de patrón `-atan2(Im(Cm),Re(Cm))/m`. No se promedian amplitudes o fases locales. La fase pierde significado cuando el modo está al nivel del redondeo; los tests de fase exigen amplitud significativa. La máscara vacía produce coeficientes/amplitudes/fases cero. El tiempo se conserva en su columna.

Phi permanece completo por proceso según T01. Con nz=1 se omite la serie PPI, igual que en la ruta serial. El diagnóstico mantiene la convención de fase existente; el desenvolvimiento temporal de fases y Fourier siguen siendo análisis de las series, sin cambios en este incremento.

## Proxy cuadrupolar

Cada bloque reutiliza `calc_gw_strain(...,local_stress)` para integrar las tres componentes Finn–Evans existentes: Ixx_ddot, Iyy_ddot e Ixy_ddot. Las coordenadas, la transformación de velocidades eulerianas a tasas coordenadas/cartesianas, el Jacobiano plano KS y su factor radial logarítmico son los mismos de la ruta serial. EF utiliza spin geométrico cero. Se conservan la resta de atmósfera y la exclusión de materia interior al horizonte.

Después de sumar esas componentes hacia rango 0:

    h_plus  = (Ixx_ddot - Iyy_ddot)/distance,
    h_cross = 2*Ixy_ddot/distance,
    distance = 1000 (unidades geométricas del código).

No se introducen derivadas temporales numéricas nuevas. Es el mismo proxy de campo débil, con proyección polar y normalización anteriores; no es extracción de ondas de una métrica dinámica. La interpretación de cortes sagitales conserva las limitaciones del diagnóstico serial y no equivale a una señal física 3D validada.

Para señales axisimétricas casi nulas, el error absoluto de h se normaliza con la escala de las componentes del tensor antes de proyectar, dividida por la distancia. Dividir por la señal cancelada mediría error relativo respecto de cero. Se conserva el umbral de 1e-11 y se comparan también las componentes sin proyectar. El fixture publica `reference_integrals.dat` para que esa escala sea auditable.

## Acreción

La definición conservada es

    mdot = -sum_jk D(i_probe,j,k)*(alpha*v^r-beta^r)*dy*dz.

`idx_probe_global` identifica el radio y cada rango verifica su `idx_probe` local: índice válido si contiene ese radio y cero si no lo contiene. Los bloques de la columna radial propietaria integran sus respectivos segmentos polares; los demás aportan cero. La suma hacia rango 0 combina la esfera angular una sola vez. En radio logarítmico se mantienen la densificación y velocidad coordenada propias de la malla; no se añade otra conversión radial a esta fórmula.

El signo de entrada es positivo. Los controles prescritos incluyen flujo radial conocido y una sonda en el cuerpo, además de atmósfera; también se compara la sonda cercana al horizonte de los toros evolucionados.

## API y salida

Inicializar malla/caché/estado, runtime y `initialize_parallel_diagnostics(layout)` después de preparar los contadores. Esto mantiene los requisitos D01: CPU, `numerical_diagnostics=true`, sin restart distribuido. Llamar entre pasos RK completos:

```fortran
call write_distributed_diagnostics(integration_time)
call write_distributed_observables(integration_time)
```

Ambas llamadas son colectivas. La segunda comprueba época, controles comunes, estado físico finito, geometría y propiedad de sonda antes de elegir qué series emitir. Respeta `do_ppi_diagnostics`, `do_gw_extraction` y `do_mdot_extraction`. En Minkowski y para bh_mass<0.1, no emite GW/mdot, como el código serial.

Bajo `output_folder/scheme_name/`, sólo rango 0 abre:

| Archivo | Columnas |
| --- | --- |
| `mpi_ppi_modes.dat` | t, C0, y por m=1..4: ReCm, ImCm, amplitud normalizada, fase de patrón; 18 columnas |
| `mpi_GW_signal.dat` | t, h_plus, h_cross; 3 columnas |
| `mpi_m_dot.dat` | t, mdot; 2 columnas |

Los nombres con prefijo evitan mezclar estas series con las anteriores. Cada serie exige archivo nuevo en su primera muestra y tiempo estrictamente creciente después. No crea carpetas; el llamador debe prepararlas antes. Se propagan colectivamente errores de apertura/escritura/cierre, muestras duplicadas y archivos existentes. No hay aún continuidad transaccional entre series y checkpoints; corresponde a C01/C03.

Para reutilizar integrales sin E/S se añaden argumentos opcionales:

| Rutina | Resultado local |
| --- | --- |
| `calc_ppi_modes(t,local_coefficients)` | vector `[C0, ReC1..ReC4, ImC1..ImC4]` |
| `calc_gw_strain(t,local_stress)` | `[Ixx_ddot,Iyy_ddot,Ixy_ddot]` |
| `calc_m_dot(t,local_mdot)` | contribución local con signo de acreción; cero en no propietarios |
| `reduce_observable_sum(values)` | suma colectiva con resultado válido sólo en rango 0 |
| `ppi_mode_values(t,coefficients,values)` | normalización/fase sobre coeficientes ya globales |

Estos ayudantes no reemplazan la validación del publicador colectivo; la salida normal distribuida debe entrar por él.

## Ajuste de D01

`collect_global_state(state,root_only=.true.)` y `collect_conservation_balance(balance,root_only=.true.)` finalizan el resultado únicamente en rango 0. En otros rangos queda el tipo inicializado a cero, sin significado global. El argumento omitido mantiene los resultados globales en todos los procesos para las pruebas y llamadores anteriores. Se comprueba que todos elijan el mismo destino antes de cambiar de colectiva.

`write_distributed_diagnostics` usa ahora ese destino raíz. Sumas, mínimos, contadores y suma de cuadrados L2 se reducen hacia 0; los máximos que sirven de escala para los cuadrados L2 siguen con Allreduce, pues cada bloque los necesita antes de sumar. Las normas de error de referencia siguen disponibles en todos los rangos. No se alteran los acumuladores ni pesos RK de D01.

## Verificación y reproducción

Cierre: 960estados prescritos y128comparaciones de evolución en los dos builds MPI; cuarenta rechazos coordinados. Campos y dt exactos. Errores máximos analítico2.172e-14, componentes8.997e-15 y series1.695e-12, todos por debajo de1e-11. 32regresiones Fortran ydos grupos de integración por cada uno de cuatro builds sinMPI/MPI opt/debug; 544comparaciones/60rechazos D01 adicionales conservan sus criterios de masa/J. 28controles del programa habitual reproducen porbytes campos, contadores y presupuesto frente a pre-D02 con iguales opciones de arquitectura. Las series conservan formatos/cadencia y equivalencia numérica. PPI optimizado cambia últimos bits (máximo error normalizado1.404e-16) por la especialización LTO (sin/cos vectorial anterior frente a sincos escalar actual), comprobada en ensamblados; se registra su error con el criterio1e-11 ya fijado. En debug PPI y todas las demás series siguen idénticas porbytes. Las pruebas de checkpoint/salidas con rutas temporales históricas fijas se repiten secuencialmente en los cuatro builds. Hashes verifican70fuentes actuales,52archivos protegidos del padre y41referencias inmutables.

Los controles usan FP64, dos hilos por rango y grids desiguales 1×1, 2×1, 1×2, 4×1, 1×4, 2×2, 4×2 y 2×4 sobre un solo equipo local.

- Prescritos: malla37×19×32, KS a=0/0.9 con radio lineal/log y EF log; doce muestras por configuración. Nueve estados con modos m1/m2 rotantes y crecimiento conocido, y controles axisimétrico, atmósfera y materia interior al horizonte. Coeficientes y combinación TT contrastados con integrales azimutales analíticas. Ghosts NaN excluidos. Sonda dentro del cuerpo y norma L2/contadores comparados entre destinos raíz/todos.
- Evoluciones: FD/FM a=0/0.9 sagitales y 3D, mallas37×19×1/8 hasta0.1M, MP5 con a=0 y WENO5 con a=0.9, sensor activo; HLLE para FD y HLLC para FM. Ocho controles por grid/build, estados iniciales copiados por índice desde referencia completa; campos y secuencia dt deben ser exactos. No se cierra D03 ni equilibrio prolongado.
- Publicación: una muestra por tiempo y archivo, esquema/normalización/fase/signos frente a salida serial. Diez fallos deliberados con2/4procesos: controles/espín, propiedad de sonda, estado no finito, tiempo equivocado, RK parcial, destino raíz discordante, duplicado, archivo existente y E/S fallida. No se aceptan timeouts como rechazos válidos.
- Regresión: suite completa sinMPI/MPI optimizada/debug con bounds/traps; presupuesto D01 multiproceso y equivalencia del programa habitual frente a binarios pre-D02.

```bash
make test-observables
make USE_MPI=1 test-observables-mpi
make USE_MPI=1 test-diagnostics-mpi
make test
```

Ejecutar matrices MPI una a la vez, sin `-j` entre las matrices. El harness usa carpetas únicas y `--bind-to none`, con espera OpenMP pasiva. Las comparaciones no miden escalabilidad.

Los primeros intentos del fixture imponían flujo casi estacionario dentro del horizonte/ergorregión y fallaron su prueba causal; se corrigió sólo ese dato prescrito a reposo euleriano para r<=3. Se preservan los logs. También se conserva la primera comparación de h casi cero y la revisión de su escala con tensor explícito. No se cambió el solver para acomodar estos controles.

Siguiente tras verificar D02: D03, perturbación reproducible por índices globales e inyección única. Checkpoints/observador/E/S MPI de producción, continuidad real24h, GPU y escalabilidad multinodo siguen pendientes. No hay nuevas consultas/cambios a CADS ni copia a D, Git, TAMU o desarrollo magnético en este incremento.
