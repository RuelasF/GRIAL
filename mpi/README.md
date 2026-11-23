# GRIAL hidrodinámico con MPI y OpenMP

Esta carpeta contiene la implementación hidrodinámica de referencia de GRIAL en `main`.
Compila con OpenMP para un proceso o con MPI y OpenMP para varios procesos.
Para toros 3D,
`--grid NR NTHETA` reparte radio y ángulo polar y conserva el azimut completo en
cada bloque. El producto de los dos números debe coincidir con el número de rangos.
Los casos 1D y el jet x–z se reparten en el primer eje.

## Descargar, compilar y probar en Linux

Instala GNU Fortran, Make, Python 3 y Open MPI con su paquete de desarrollo
(`gfortran make python3 openmpi-bin libopenmpi-dev` en Debian y Ubuntu).
Desde una cuenta de usuario normal, ejecuta lo siguiente.

```bash
git clone --single-branch --branch main https://github.com/RuelasF/GRIAL.git
cd GRIAL/mpi
make -j2 USE_MPI=1 ARCH_FLAGS=''
OMP_NUM_THREADS=1 mpirun -n 2 ./grial_mpi --distributed benchmarks/mpi_reference/sod.par --grid 2 1
```

El caso Sod escribe diagnósticos, piezas VTK y un checkpoint en
`native/weno5_hlle/`.
Para otra corrida, selecciona un nuevo `output_folder` en el archivo de parámetros.
`OMP_NUM_THREADS` controla los hilos **por rango**. Esta entrada requiere
`numerical_diagnostics=true`, cadencias de diagnósticos y extracción positivas
y que el presupuesto radial de la entrada serial esté desactivado.

Para compilar y ejecutar en un solo proceso, bastan GNU Fortran y Make.
Usa un directorio de salida distinto si ya ejecutaste el caso MPI.

```bash
make -j2 ARCH_FLAGS=''
OMP_NUM_THREADS=4 ./grial benchmarks/mpi_reference/sod.par
```

Para probar compilación y ejecución con Sod, toro 3D y reinicio con otro reparto:

```bash
OMP_NUM_THREADS=1 make USE_MPI=1 ARCH_FLAGS='' test-application
```

Consulta la [guía de la entrada distribuida](docs/ENTRADA_DISTRIBUIDA_CONTROL_C.md)
para validar parámetros con `--check` y continuar con `--resume`. La
[guía de contenedores](containers/README.md) describe el Dockerfile CPU.
La [procedencia](PROVENANCE.json) registra las fuentes hidrodinámicas de partida.

## Elegir el intercambio MPI

Esta versión incorpora las fuentes de la copia hidrodinámica optimizada. El
intercambio bloqueante sigue siendo la configuración de referencia.

| Compilación | Ejecutable | Intercambio |
| --- | --- | --- |
| `make USE_MPI=1 OVERLAP=0 ARCH_FLAGS=''` | `./grial_mpi` | Bloqueante, con `MPI_Sendrecv` |
| `make USE_MPI=1 OVERLAP=1 ARCH_FLAGS=''` | `./grial_mpi_overlap` | No bloqueante, con solapamiento experimental |

Ambas compilaciones pueden coexistir. Lanza el ejecutable elegido con los
mismos argumentos MPI, archivo de parámetros y reparto. Usa un directorio
de salida nuevo para cada corrida. La variante con solapamiento adelanta
cálculo interior del jet y los toros FD/FM 3D mientras intercambia halos
radiales y, cuando hay reparto polar, también polares. En CADS el efecto en
rendimiento varió según el problema y el reparto.

También puedes añadir `CFL_CACHE=1` para reutilizar el límite causal del paso
de tiempo calculado sobre la geometría fija. La opción viene desactivada y
agrega el sufijo `_cfl_cache` al ejecutable. `PROFILE=1` construye el perfilador
`grial_mpi_profile` y requiere `OVERLAP=0`.

Con el perfilador activo, `GRIAL_COLLECTIVE_TRACE=1` escribe
`grial_collectives_step_<inicio>_rank_<rango>.tsv` junto a los diagnósticos de
la corrida. El nombre usa el paso inicial de la sesión. Si la traza ya existe
o falla su escritura, el perfilador avisa y la simulación puede publicar su
checkpoint final. Cada rango registra
hasta 131072 llamadas con el sitio y las marcas de entrada y retorno.
En un mismo host GNU/Linux se pueden comparar esas marcas entre rangos.
El intervalo posterior a la última entrada también incluye el trabajo local
del colectivo. Las marcas de distintos nodos carecen de una referencia común.

### Malla polar para toros 3D

En `--distributed`, los toros FD/FM en Kerr–Schild admiten
`polar_mesh_h = 0.3` en el archivo de parámetros. Para la coordenada lógica
uniforme `q`, el ángulo físico es `theta = q + (1-h) sin(2q)/2`.
El valor predeterminado `h = 1` conserva la malla uniforme. Con `0 < h < 1`
hay más resolución angular cerca del ecuador. El checkpoint guarda el mapa
para `--resume` y las piezas VTK se escriben en coordenadas físicas.
Un checkpoint mapeado requiere un ejecutable que admita esta opción.

Para actualizar la versión de referencia, apunta tu clon a `main`, ejecuta
`git pull --ff-only` y recompila la variante deseada desde `mpi/`.

La regresión `make USE_MPI=1 ARCH_FLAGS='' test` puede ejecutarse en el clon
público. Las autopruebas del controlador de campañas largas de SRHD requieren
el script de Slurm de la instalación local.

## Checkpoints distribuidos

C01 guarda bloques físicos UP/P, configuración global, contadores y época; publica un manifiesto y un marcador final después de comprobar tamaños y leer los datos guardados. Una escritura incompleta deja la generación sin publicar y conserva la anterior. La recuperación al mismo reparto comprueba todos los bloques antes de modificar la memoria. [API, formato, pruebas y límites](docs/CHECKPOINT_DISTRIBUIDO_C01.md).

```bash
make USE_MPI=1 test-checkpoint-mpi
```

C02 carga configuración antes de reservar la malla, lee intersecciones con otro reparto e importa serial V7/V8/V9 sin reunir el campo global. [API, convención de contadores y pruebas](docs/REINICIO_REPARTO_C02.md).

```bash
make USE_MPI=1 test-restart-mpi
```

C03 restaura presupuestos/observadores y conecta segmentos mediante un driver separado. Guarda un sidecar y marcador final ligados a C01, valida las historias antes de recuperar y conserva el calendario sin duplicar muestras. [Contrato, comandos y resultados](docs/OBSERVADORES_SEGMENTOS_C03.md). El V9 serial conserva todos sus registros.

```bash
make USE_MPI=1 mpi-segment-driver
make USE_MPI=1 test-observer-mpi
```

448 comparaciones segmentadas/continuas con campos exactos y error de series ≤5.499e−15; cuatro builds optimizados/debug sin/con MPI aprobados. Control C está integrado localmente; validación integral y escalabilidad siguen pendientes.

## Descomposición de dominio

[Reparto, variables globales/locales, pruebas y límites M02](docs/DESCOMPOSICION_M02.md). Se parte de SRHD1D, después SRHD2D y finalmente GRHD. M03 intercambia halos con propietarios físicos, M04 verifica reconstrucción/RHS, M05 coordina CFL/eventos y M06 evoluciona bloques con RK mediante drivers. T01 verifica el reparto r×theta, con phi completo, mediante ocho grids de 1/2/4/8 procesos y controles opt/debug. [Contrato, comandos y límites](docs/DOMINIO_RTHETA_T01.md). D01 verifica balances y contadores globales; [API, unidades y series](docs/DIAGNOSTICOS_GLOBALES_D01.md). D02 verifica PPI/acreción/ondas con integración local, reducciones hacia rango0 y finalización/salida allí; [API y pruebas](docs/OBSERVABLES_D02.md). D03 verifica umbral global, secuencia histórica de ruido, evento único y presupuesto explícito; [API y límites](docs/PERTURBACIONES_D03.md). C01 guarda; C02 carga configuración y recupera con otro reparto/importa serial; C03 continúa presupuestos/series entre sesiones. Siguiente V01.

```bash
make -j2 USE_MPI=1 ARCH_FLAGS='' test-domain-mpi
```

## Intercambio de primitivas

[API, fronteras y evidencia M03](docs/HALOS_M03.md). Buffers persistentes, dos Sendrecv por eje radial y llenado de aristas/esquinas tras recibir. Comparado con fronteras seriales y fuentes anteriores, con error0. El integrador distribuido se usa mediante drivers CPU y la entrada explícita de Control C.

```bash
make -j2 USE_MPI=1 ARCH_FLAGS='' test-halos-mpi
```

## Reconstrucción e interfaces

[Metadatos globales, criterios y límites M04](docs/INTERFACES_M04.md). Se preserva CPU serial, se alinean las políticas por cara y se comprueban estados/flujos/RHS con1/2/4procesos. Control C habilita su entrada explícita, tras verificar la integración local.

```bash
make -j2 USE_MPI=1 ARCH_FLAGS='' test-interfaces-mpi
```

## CFL y eventos

[Fórmula, API, cadencia y evidencia M05](docs/CFL_EVENTOS_M05.md). Primero se obtienen los tres máximos globales por dirección, después se forma dt. Se verifica época/calendario y se propagan solicitudes de parada/checkpoint. La política serial se conserva; estos controles no habilitan todavía evolución del fluido con varios procesos.

```bash
make -j2 USE_MPI=1 ARCH_FLAGS='' test-time-mpi
```

## Evolución Runge–Kutta

[Rutina compartida, controles y límites M06](docs/RK_DISTRIBUIDO_M06.md). Cien trayectorias con cinco reconstructores y HLLE/HLLC, desde ondas/choques 1D hasta pilotos Michel y FD/FM, con 1/2/4 procesos. Campos y secuencia dt exactos; se conserva el programa serial frente a M05. Los drivers no utilizan los escritores de producción.

```bash
make -j2 USE_MPI=1 ARCH_FLAGS='' test-rk-mpi
```

## Choques hasta t=0.4

[Perfiles, diferencias y referencia exacta M06B](docs/CHOQUES_MPI_T04_20261024.md). Sod y choque fuerte con N=400/401: 128 trayectorias, campos y secuencias temporales exactos con MPI 1/2/4, incluidos controles debug. Se presentan por separado los errores de discretización y el diagnóstico de borde de Godunov.

```bash
make -j2 ARCH_FLAGS='' analysis/mpi_shock_profiles
make -j2 USE_MPI=1 ARCH_FLAGS='' analysis/mpi_shock_profiles
python3 tools/compare_shocks_mpi.py --help
```

Este driver de análisis escribe archivos propios por rango. No habilita las salidas ni reinicios MPI del programa principal.

## Shu–Osher con MPI 1/2/4/8

[Perfiles, ampliaciones y diferencias M06C](docs/SHU_MPI_20261024.md). 80 trayectorias hasta t=0.4, con N=400/401: datos iniciales, campos finales y trazas idénticos con 1/2/4/8 procesos, incluidos controles debug. Se comparan cinco reconstructores y HLLE/HLLC usando el preset de GRIAL y el mismo núcleo de M06.

```bash
python3 tools/compare_shu_mpi.py --help
```

## Compilar la base actual

```bash
make -j2 ARCH_FLAGS=''
make test ARCH_FLAGS=''
./grial --check benchmarks/mpi_reference/fd_3d_a09.par
```

Los archivos de `benchmarks/mpi_reference/` son controles cortos locales, no una campaña de producción ni scripts de envío al cluster. Para ejecutarlos, emplear una carpeta de resultados nueva.

La referencia serial se conserva en `reference/grial_hydro_base_20261020.tar.gz`, con SHA256 y manifiesto. El archivo comprimido contiene las fuentes de referencia. El verificador también comprueba los parámetros de `reference/controls/`. Para comprobar la referencia:

```bash
python3 tools/check_reference.py
```

El modo `--working-copy` compara los archivos de trabajo con la base inicial. Es apropiado al cerrar MPI00; mostrará diferencias esperables después de comenzar modificaciones MPI. No restaura ni sobrescribe archivos.

## Procedencia y límites

Se copian fuentes hidrodinámicas congeladas de la campaña FD/FM F03, anteriores a conectar SRMHD al ejecutable. Sus módulos de condiciones, ecuaciones, flujos, evolución, métricas, salida y aceleración coinciden porbytes con la versión principal al preparar esta copia. Las diferencias de interfaz se documentan en la evidencia de creación. Se conservan las23pruebas hidrodinámicas; no se copian módulos ni pruebas SRMHD.

La base heredada contiene una bandera interna reservada `is_mhd=false`; no hay aquí una ruta magnética operativa ni se admite activar SRMHD desde los parámetros. Los campos evolucionados siguen siendo las cinco variables hidrodinámicas.

Los enlaces a expedientes de validación locales en las notas históricas indican la procedencia de las pruebas.

## Observables distribuidos

[D02, formatos y alcance](docs/OBSERVABLES_D02.md). Las integrales recorren los bloques; rango0 normaliza PPI, proyecta el tensor y publica las tres series. D01 permite destino raíz para salida y conserva los colectores globales por defecto. No se reúne la malla. Control C conecta esta API a la entrada MPI explícita.

```bash
make USE_MPI=1 test-observables-mpi
```


## Visualización distribuida

C04 escribe una pieza XML `.vtu` por rango, un índice `.pvtu` y el tiempo `.pvd`, sin reunir la malla. Conserva las caras físicas y prismas en los polos. El driver de segmentos escribe inicial/periódica/terminal y continúa sin repetir el frame de reinicio. [Formato, lector, colección y resultados](docs/VISUALIZACION_BLOQUES_C04.md).

```bash
make USE_MPI=1 test-visualization-mpi
python tools/parallel_vtk.py data/weno5_hlle/mpi_fields --collection data/weno5_hlle/mpi_fields/frames.pvd
```

Abrir `frames.pvd` en ParaView para recorrer el tiempo. El lector/inspector usa sólo Python estándar; reconstruir campos globales es opcional y posterior. Verificado con VTK9.7.1 y234controles de geometría,40cadenas de evolución y suites opt/debug. No certifica GUI, multinodo o rendimiento de miles de piezas; Control C conecta esta salida a la entrada MPI explícita; V01–V04 siguen pendientes.
