# D03: perturbación global reproducible y evento único

Implementación local cerrada el 25/oct/2026 en la copia hidrodinámica `mpi/`. [Resultados verificados](../.validation/d03_perturbations_20261025_v1/VERIFIED.json); el expediente conserva parámetros, logs, binarios, fuentes anteriores y parche.

## Contrato y selección de materia

`parallel_perturbations:apply_parallel_perturbation(layout,injected)` es colectivo. Todos los procesos deben llamarlo en el mismo punto del calendario, fuera de una etapa RK, con el dominio, halos y diagnósticos D01 preparados. El argumento lógico opcional `injected` indica si esa llamada aplicó el evento. El estado físico debe ser finito y causal, y los controles, tiempo, paso y bandera de inyección deben coincidir entre rangos.

El umbral conserva la definición histórica:

\[
\rho_{\max}=\max_{\text{celdas físicas globales}}\rho,
\qquad \mathcal B=\{(i,j,k):\rho_{ijk}>0.05\rho_{\max}\}.
\]

Se excluyen los ghosts; una celda exactamente en el 5% queda fuera. Un rango puede no contener ninguna celda de la selección. Usar máximos locales produciría máscaras y perturbaciones distintas entre repartos, por lo que D03 utiliza el máximo global de D01.

## Compatibilidad de las tres perturbaciones

Se reutilizan las fórmulas de `conditions.f90`, ahora con argumentos opcionales para el máximo global y, en los dos tipos de ruido, un campo de uniformes proporcionado por el llamador:

- Presión: \(p'=p[1+A(2u-1)]\), sin cambiar densidad o velocidad.
- Densidad aleatoria: \(f=1+A(2u-1)\), \(\rho'=f\rho\), \(p'=f^\Gamma p\); conserva la relación isentrópica local.
- Modo de densidad: \(f=1+A\cos(m\phi)\), con la misma actualización isentrópica y las coordenadas globales existentes.

Las tres seleccionan la materia mediante la densidad anterior a la inyección. Los valores fuera del cuerpo y las velocidades se conservan. Se admiten amplitudes \(0\leq A<1\). Una amplitud cero sigue siendo un evento y consume la bandera de inyección, igual que la ruta habitual.

La llamada histórica de tres argumentos sigue disponible para los ruidos. Al omitir las entradas nuevas, las rutinas mantienen la selección, semilla y generación escalar anteriores. `set_reproducible_random_seed` forma el vector de semillas con enteros de 64 bits antes de reducirlo al rango entero del generador; evita el desbordamiento de semillas extremas. Las semillas normales conservan la secuencia previa. El resultado depende de la implementación intrínseca del compilador: no se promete identidad del generador entre compiladores o versiones de biblioteca distintos.

## Orden global del ruido

El generador intrínseco histórico consume un uniforme por celda seleccionada, recorriendo `k`, después `j` y finalmente `i`. D03 conserva esa secuencia:

1. Cada rango construye su máscara física con el máximo global.
2. Rango 0 reúne las máscaras y las coloca mediante los offsets globales del reparto desigual `x/r × y/theta`.
3. Rango 0 inicializa la semilla y recorre la máscara completa en el orden histórico. Las celdas no seleccionadas no consumen números.
4. Empaqueta los uniformes por bloque y los distribuye; cada rango aplica la misma rutina de `conditions` a sus primitivas.

No se reúne el estado del fluido. En un solo proceso se omiten los gather/scatter y se conserva el mismo recorrido. La tercera dirección permanece completa, conforme a T01. El modo coherente no requiere generar ni distribuir números aleatorios.

Esta primera solución prioriza conservar los datos iniciales históricos. La memoria temporal de rango 0 crece con el número global de celdas; el gather/scatter usa conteos MPI clásicos y rechaza tamaños superiores al máximo entero admitido. No constituye una medición de escalabilidad ni un método de generación distribuida de memoria constante.

## Calendario y protección contra repetición

Se reutiliza `coordinate_time_actions` de M05. La perturbación se aplica al inicio del primer paso cuyo tiempo ya alcanzó `perturbation_time`; no se modifica `dt` para caer exactamente sobre el umbral. Si la simulación ya alcanzó `final_time`, prevalece el cierre y no se inyecta. El driver calcula `dt` antes de la perturbación, en el mismo orden que el programa habitual.

Tras cambiar las primitivas se recalculan `U` y `UP` en las celdas físicas y se rellenan los halos. Todos los rangos establecen `perturbation_applied=.true.`. Las siguientes llamadas son no operativas. Desacuerdos de bandera, semilla o reloj se rechazan colectivamente.

El formato serial V9 ya conserva esa bandera. Las pruebas D03 escriben y leen V9 antes y después del evento para verificar su efecto al repartir un estado controlado. **No prueban un reinicio distribuido completo**: el guardado, redistribución y restauración del presupuesto/observadores MPI corresponden a C01–C03. D01 continúa rechazando `do_restart=true` y la entrada multiproceso de producción de `grial` conserva su guarda.

## Presupuesto y recibo de inyección

La modificación deliberada de las conservadas es un término propio, separado de fuentes geométricas y de reparaciones de primitivas:

\[
\Delta Q_{\mathrm{pert}}
=\sum_{\text{celdas físicas}}(U^+ - U^-)\,\Delta x\,\Delta y\,\Delta z,
\]

\[
R=Q-Q_0-\int(F_{\mathrm{inferior}}-F_{\mathrm{superior}})\,dt
-\int S\,dt-\Delta Q_{\mathrm{rec}}-\Delta Q_{\mathrm{pert}}.
\]

`U` ya está densitizada; el volumen anterior es el de las coordenadas lógicas. El cambio se mide respecto de `UP` antes de reconstruir `U/UP`, incluyendo el pequeño ajuste de consistencia de la conversión. La inyección es instantánea y no recibe pesos RK. Los pesos anteriores del RHS y de las reparaciones no cambian.

`conservation_balance%perturbation(5)` y `accumulate_perturbation_budget` implementan ese término. `mpi_conservation_budget.dat` pasa de 57 a **62 columnas**; el orden es tiempo/paso, Q inicial (5), Q actual (5), flujo (30), fuentes (5), reparaciones (5), perturbación (5) y residuo (5). El bloque de perturbación ocupa las columnas 53–57 contando desde 1; el residuo ocupa 58–62. Los lectores deben usar el esquema nuevo o localizar el residuo como las cinco columnas finales.

Sólo rango 0 crea `output_folder/scheme_name/mpi_perturbation_event.dat`, con una fila y 13 valores: tiempo, paso, tipo, semilla, amplitud, modo, máximo global, número de celdas seleccionadas y cinco cambios integrales en orden `D,tau,S1,S2,S3`. La apertura exige un archivo nuevo; un recibo existente o una falla de E/S provoca rechazo colectivo y conserva el archivo previo. Las series seriales anteriores no cambian.

## Controles y reproducción

El programa `tests/test_parallel_perturbations.f90` compara la malla completa y bloques desiguales, verifica repetición, saldos y observables. Incluye dos campos sintéticos, materia localizada en un subconjunto de rangos, igualdad exacta en el umbral, tres perturbaciones, semillas 3435/17/extrema, amplitud cero, FD/FM con espines 0/0.9 y cortes sagitales/3D. Las evoluciones de toros alcanzan 0.04M con umbral 0.012345M; son controles del reparto y del evento, no estudios de estabilidad o PPI tardía. Las pruebas V9 controladas completan 54 casos por grid.

`tests/test_perturbations_mpi.py` ejecuta los grids `1×1`, `2×1`, `1×2`, `4×1`, `1×4`, `2×2`, `4×2` y `2×4`, y once errores deliberados con 2/4 rangos. Comprueba archivos, rechazo sin bloqueo y preservación del marcador de E/S existente. Cada ejecución necesita un directorio nuevo.

```bash
make test-perturbations
make USE_MPI=1 test-perturbations-mpi
```

El expediente contiene cuatro builds aislados, optimizado/debug sin MPI/con MPI; debug usa `-fcheck=all -ffpe-trap=invalid,zero,overflow`. `verify_increment.py` conserva las suites completas y matrices D01/D02/D03; `compare_native.py` compara el programa habitual con binarios anteriores de iguales flags. Parámetros, logs, fuentes anteriores, parche y hashes quedan junto a los resultados. Los resultados finales y sus máximos se registran en `VERIFIED.json` al cerrar el incremento.

## Resultado del cierre

Con GNU Fortran 14.2.0 y OpenMPI 5.0.7 locales, ocho grids por build optimizado/debug: **864 comparaciones** de 54 casos por grid, todas con campos idénticos. Los máximos normalizados fueron:

| Cantidad | Error máximo | Criterio |
| --- | ---: | ---: |
| Campos | 0 | identidad |
| Residuo del presupuesto con perturbación | 1.2251e−14 | <1e−11 |
| Término de inyección | 4.2964e−17 | <1e−11 |
| Observables después del evento | 4.9889e−13 | <1e−11 |

Los 44 errores deliberados se rechazaron sin bloqueo; se comprobaron 51 recibos por ejecución positiva y tres casos V9 ya aplicados sin un nuevo recibo. Las cuatro suites sin MPI/con MPI, optimizado/debug, pasaron sus 33 programas Fortran, dos grupos de integración y ocho autopruebas del controlador. Se repitieron las matrices completas D01/D02: 544 comparaciones de presupuesto y 960 estados prescritos más 128 evoluciones de observables. Su error máximo de series fue 1.6281e−12, por debajo de 1e−11.

En **120 comparaciones del programa habitual** contra binarios pre-D03 de iguales flags e hilos, las cinco conservadas, tiempo/paso, contadores, balance radial y series coincidieron por bytes. Incluyen FD/FM con a=0/0.9, sagital/3D y las tres perturbaciones, además de seis controles sin inyección. La referencia completa de 41 archivos y los 193 archivos protegidos, incluidos fuentes/binario del padre y los dos paquetes CADS, permanecieron intactos. El fallo inicial del fixture por escribir `1e-200` en precisión simple se conserva; se corrigió el literal a FP64 sin modificar tolerancias ni fórmulas físicas.

Graphify reutilizó el grafo existente de 906 nodos para ubicar las rutinas históricas; la consulta tuvo un presupuesto de salida de 1500 tokens, sin reconstrucción ni inferencia LLM adicional. El grafo no describe los incrementos MPI recientes: sus contratos se comprobaron leyendo los fuentes y ejecutando estos controles.

Las fuentes de las campañas CADS KHI/Jets v1/v2 permanecen congeladas y no reciben D03. El siguiente incremento es C01, una generación de checkpoint distribuido completa y recuperable.
