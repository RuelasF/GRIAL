ARG BASE_IMAGE=debian:bookworm-slim
FROM ${BASE_IMAGE} AS dependencies

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       libgfortran5 libgomp1 openmpi-bin python3 \
    && rm -rf /var/lib/apt/lists/* \
    && useradd --uid 1000 --create-home grial \
    && install -d -o 1000 -g 1000 /opt/grial /work \
    && chmod 1777 /work

ENV OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=PASSIVE \
    PYTHONDONTWRITEBYTECODE=1 HOME=/tmp

FROM dependencies AS build
ARG BUILD_JOBS=2
RUN apt-get update \
    && apt-get install -y --no-install-recommends gfortran make libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/grial
COPY --chown=1000:1000 Makefile *.f90 ./
COPY --chown=1000:1000 tests/*.f90 tests/*.py ./tests/
COPY --chown=1000:1000 tools/*.f90 tools/*.py ./tools/
COPY --chown=1000:1000 containers/smoke.py containers/sod.par ./containers/
USER 1000:1000

# Empty ARCH_FLAGS avoids inheriting the build computer's instruction set.
RUN make -j"${BUILD_JOBS}" ARCH_FLAGS='' all \
    && make -j"${BUILD_JOBS}" USE_MPI=1 ARCH_FLAGS='' all \
    && python3 containers/smoke.py --serial /opt/grial/grial \
       --mpi /opt/grial/.build/mpi/grial_mpi --output /opt/grial/image_smoke \
    && sha256sum Makefile *.f90 tests/*.f90 tests/*.py tools/*.f90 tools/*.py \
       containers/smoke.py containers/sod.par > SOURCE_SHA256.txt \
    && sha256sum grial .build/mpi/grial_mpi > BINARY_SHA256.txt \
    && dpkg-query -W > PACKAGES_BUILD.txt \
    && gfortran --version > TOOLCHAIN.txt \
    && mpirun --version >> TOOLCHAIN.txt

FROM dependencies AS runtime
COPY --from=build /opt/grial/grial /opt/grial/.build/mpi/grial_mpi /usr/local/bin/
COPY --from=build /opt/grial/SOURCE_SHA256.txt /opt/grial/BINARY_SHA256.txt \
     /opt/grial/PACKAGES_BUILD.txt /opt/grial/TOOLCHAIN.txt /opt/grial/
COPY --from=build /opt/grial/image_smoke/SUMMARY.json /opt/grial/SMOKE_BUILD.json
COPY tests/*.py /opt/grial/tests/
COPY tools/*.py /opt/grial/tools/
COPY containers/smoke.py containers/sod.par /opt/grial/containers/
RUN dpkg-query -W > /opt/grial/PACKAGES_RUNTIME.txt
USER 1000:1000
WORKDIR /work
CMD ["grial_mpi", "--parallel-check"]
