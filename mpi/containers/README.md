# GRIAL CPU: trasladar un entorno OpenMP/MPI con Docker

Este paquete corresponde a la copia hidrodinámica `mpi/`. Incluye `grial` y `grial_mpi`, sus bibliotecas GNU/OpenMP/OpenMPI y Python estándar para lectores/pruebas. Compila con `ARCH_FLAGS=''`, FP64 y los flags habituales del Makefile. No contiene la versión magnética, GPU, resultados ni campañas congeladas de CADS.

El receptor necesita Docker funcionando y una arquitectura compatible. Una imagen `linux/amd64` sirve para equipos Linux x86-64 Intel/AMD compatibles con esa base; no se convierte en ARM por quitar `-march=native`. No requiere instalar Fortran, MPI o Python en el host. La construcción inicial necesita Internet; trasladar el TAR permite cargar esa misma imagen sin descargar dependencias.

**Estado27/oct/2026:** Docker/Podman/Apptainer no están instalados en esta sesión. Compilación portable local sin/conMPI con GNU14.2/OpenMPI5.0.7 y prueba de cuatro perfiles aprobadas: campos exactos y visualización válida. Es una ejecución nativa sobre las fuentes seleccionadas, no dentro de Debian bookworm/Docker. **Imagen Docker pendiente de construir y probar**, al igual que el equipo receptor. No se instaló un daemon ni se modificaron CADS/jobs. V01 permanece independiente.

## 1. Construir en una computadora con Docker

Desde la copia actual del proyecto:

```bash
cd ~/Codigo/GRIAL/mpi
docker build --build-arg BUILD_JOBS=2 --tag grial-mpi:cpu-20261027 .
docker run --rm grial-mpi:cpu-20261027
docker image inspect --format '{{.Id}} {{.Os}}/{{.Architecture}}' grial-mpi:cpu-20261027
```

También se preparó un paquete de fuentes pequeño en `.validation/docker_portable_20261027_v1/grial-mpi-source_20261027.tar.gz`, con su SHA256. **Ese TAR.GZ contiene fuentes, no una imagen Docker**: permite construirla en otra computadora cuando aquí no haya Docker. Tras copiarlo, usar:

```bash
sha256sum --check grial-mpi-source_20261027.tar.gz.sha256
tar -xzf grial-mpi-source_20261027.tar.gz
cd grial-mpi-source_20261027
docker build --build-arg BUILD_JOBS=2 --tag grial-mpi:cpu-20261027 .
```

La compilación se realiza dentro de una etapa `build`, con Debian bookworm, GNU Fortran y OpenMPI de esa distribución. Se construyen los ejecutables sin/conMPI y una prueba real de Sod debe pasar antes de generar la etapa `runtime`. La imagen final lleva los ejecutables/bibliotecas y lectores; los compiladores y objetos quedan en `build`.

La prueba usa OpenMP1 y MPI1/2/4, un hilo por proceso, 64celdas hasta t=0.02, WENO5–HLLE/sensor. Comprueba parámetros sin escritura, estado final, tamañoMPI, U/P, contadores físicos y piezas/colección VTK. Umbral1e−11; los procesos pueden sobresuscribirse sólo en ese control corto. No es un benchmark ni sustituye las suites científicas o V01 multinodo.

Se guardan en `/opt/grial` los recibos `SMOKE_BUILD.json`, `TOOLCHAIN.txt`, `PACKAGES_BUILD.txt`, `PACKAGES_RUNTIME.txt`, `SOURCE_SHA256.txt` y `BINARY_SHA256.txt`. Los manifiestos de fuente/binarios usan las rutas originales de compilación. Los paquetes/base no están fijados a un snapshot de repositorio: guardar el TAR conserva el artefacto exacto, mientras que una reconstrucción posterior puede resolver versiones diferentes. `BASE_IMAGE` admite un digest mediante `--build-arg` si se necesita fijarlo.

El comando por defecto hace `grial_mpi --parallel-check` con un proceso; debe terminar mostrando MPI activo. No lanza una simulación automáticamente.

## 2. Guardar y transportar

En la computadora donde se construyó:

```bash
docker save --output grial-mpi_cpu_20261027.tar grial-mpi:cpu-20261027
sha256sum grial-mpi_cpu_20261027.tar > grial-mpi_cpu_20261027.tar.sha256
```

Copiar **TAR y SHA256** por USB/SCP a la computadora receptora. No hace falta publicar la imagen en Docker Hub. El TAR contiene las capas y dependencias; las simulaciones y sus `.par` se trasladan aparte.

En el receptor, dentro de la carpeta con ambos archivos:

```bash
sha256sum --check grial-mpi_cpu_20261027.tar.sha256
docker load --input grial-mpi_cpu_20261027.tar
docker run --rm grial-mpi:cpu-20261027
```

Si todavía no tiene Docker, usar las [instrucciones oficiales para su distribución](https://docs.docker.com/engine/install/). No ejecutar un instalador genérico sin conocer su sistema. No se necesita el código fuente en el receptor para ejecutar la imagen cargada.

## 3. Ejecutar un caso con resultados persistentes

Ejemplo Linux/bash. Usar una carpeta **nueva para cada ejecución**: GRIAL rechaza sobrescribir una corrida existente. El ejemplo exporta un `.par` incluido en la imagen; para producción sustituirlo por el caso propio.

```bash
mkdir -p ~/GRIAL_runs/sod_mpi2
cd ~/GRIAL_runs/sod_mpi2
docker run --rm grial-mpi:cpu-20261027 cat /opt/grial/containers/sod.par > caso.par
docker run --rm --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=$PWD,dst=/work" \
  -e OMP_NUM_THREADS=1 \
  grial-mpi:cpu-20261027 \
  mpirun -n 2 --map-by core --bind-to core \
  grial_mpi --distributed /work/caso.par --grid 2 1
```

Usar un usuario Linux ordinario; el ejemplo conserva su UID/GID en los archivos. El contenedor tiene usuario no-root por defecto y no desactiva la protección de OpenMPI frente a ejecución como root. Si el host usa rootless Docker/user namespaces o ACL especiales, comprobar los permisos de la carpeta compartida.

Los resultados quedan en `~/GRIAL_runs/sod_mpi2/native/weno5_hlle/`, aunque se elimine el contenedor con `--rm`. `/work` es la carpeta compartida: conservar allí parámetros, resultados y checkpoints. Para otros `.par`, usar `output_folder` relativo o dentro de `/work`; una ruta fuera del montaje no persistirá al retirar el contenedor.

Para ejecutar OpenMP en una **carpeta distinta**, exportar/copiar el `.par` y sustituir la última parte por:

```bash
docker run --rm --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=$PWD,dst=/work" \
  -e OMP_NUM_THREADS=4 \
  grial-mpi:cpu-20261027 \
  grial --distributed /work/caso.par
```

Los números 2/4 son ejemplos: respetar los cores físicos disponibles y las tareas de otros usuarios. Para híbrido, asignar suficientes cores a cada rango con `--map-by slot:PE=H --bind-to core` y `OMP_NUM_THREADS=H`, verificando `--report-bindings`. `--cpus` de Docker limita cuota; no equivale a elegir sockets/cores. Tras inspeccionar la topología puede usarse `--cpuset-cpus` para restringir CPU concretas.

## 4. Repetir el control en el receptor

En una carpeta de resultados escribible y con `smoke_receptor` aún inexistente:

```bash
docker run --rm --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=$PWD,dst=/work" \
  grial-mpi:cpu-20261027 \
  python3 /opt/grial/containers/smoke.py --output /work/smoke_receptor
```

Revisar `smoke_receptor/SUMMARY.json`: exigir `passed=true`. Los logs y campos se conservan si falla. Esto vuelve a ejecutar OpenMP/MPI y el lector sobre el equipo receptor; consultar sólo `SMOKE_BUILD.json` demostraría la construcción original, no la ejecución en el nuevo equipo.

Los reinicios usan la entrada existente `grial_mpi --distributed --resume DIRECTORIO --grid R THETA`. Tomar el directorio de `checkpoint=` en `run_status.dat` y mantener el mismo montaje `/work` y sus historias; conservar la misma imagen para una cadena de reinicio. Los controles de integridad/configuración del solver siguen activos. No prometer portabilidad binaria universal de checkpoints entre compiladores/arquitecturas distintos.

## Alcance para la computadora de dos sockets

Dos sockets son un nodo: esta imagen permite ejecutar y comparar MPI/OpenMP en ese equipo. Dos contenedores del mismo host tampoco demuestran comunicación entre máquinas físicas. El Dockerfile no prepara SSH, red, Slurm, InfiniBand o un lanzador multinodo. Eso requiere integración específica con la instalación del centro; no sustituir el job V01 ni mezclar el MPI del host con las bibliotecas de la imagen sin comprobar compatibilidad.

Para mediciones NUMA/rendimiento usar preferentemente Linux nativo con afinidad y recursos registrados. Docker Desktop utiliza una VM Linux y no ofrece por sí solo una medición equivalente del servidor físico.

Para disponer también de fuentes y compiladores dentro de una imagen de desarrollo, construir la etapa existente con otro tag: `docker build --target build -t grial-mpi:devel .`. La imagen runtime guardada arriba basta para ejecutar.

Referencias técnicas: [etapas Docker](https://docs.docker.com/build/building/multi-stage/), [save](https://docs.docker.com/reference/cli/docker/image/save/), [load](https://docs.docker.com/reference/cli/docker/image/load/), [montajes](https://docs.docker.com/engine/storage/bind-mounts/), [mapeo MPI](https://docs.open-mpi.org/en/main/man-openmpi/man1/mpirun.1.html), [Docker Desktop y VM](https://docs.docker.com/desktop/setup/install/linux/).
