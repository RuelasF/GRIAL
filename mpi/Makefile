USE_MPI ?= 0
MPIFC ?= mpifort
SRC_DIR ?= .
.DEFAULT_GOAL := all

ifneq ($(filter $(USE_MPI),0 1),$(USE_MPI))
$(error USE_MPI must be 0 or 1)
endif
ifeq ($(USE_MPI),1)
ifneq ($(filter openacc-% benchmark-jet-p100,$(MAKECMDGOALS)),)
$(error M01 supports CPU MPI only)
endif
ifneq ($(MPI_BUILD_INTERNAL),1)
MPI_FORWARD_GOALS := $(if $(MAKECMDGOALS),$(MAKECMDGOALS),all)
.PHONY: __mpi_forward $(MPI_FORWARD_GOALS)
$(MPI_FORWARD_GOALS): __mpi_forward ;
__mpi_forward:
	@mkdir -p .build/mpi/tests .build/mpi/analysis
	+$(MAKE) -C .build/mpi -f $(abspath Makefile) USE_MPI=1 MPI_BUILD_INTERNAL=1 SRC_DIR=$(CURDIR) $(MPI_FORWARD_GOALS)
	@if test -f .build/mpi/grial_mpi; then ln -sfn .build/mpi/grial_mpi grial_mpi; fi
else
MPI_BUILD_INTERNAL_BODY := 1
endif
else
MPI_BUILD_INTERNAL_BODY := 1
endif

ifeq ($(MPI_BUILD_INTERNAL_BODY),1)
# ==========================================
# Makefile para GRIAL (General Relativistic Instability Accretion Laboratory)
# ==========================================

# Compilador y Banderas de optimización
FC = gfortran
ifeq ($(USE_MPI),1)
FC = $(MPIFC)
endif
PYTHON ?= python3
ARCH_FLAGS ?= -march=native
FFLAGS = -O3 -fopenmp $(ARCH_FLAGS) -flto=auto -fno-math-errno -fno-trapping-math -ffree-line-length-none -cpp
#FFLAGS = -fopenmp -O0 -g -fbacktrace -Wall -Wextra -fcheck=all -finit-real=snan -ffpe-trap=invalid,zero,overflow,underflow,denormal

ifeq ($(USE_MPI),1)
override FFLAGS += -DUSE_MPI
ifneq ($(findstring USE_OPENACC,$(FFLAGS)),)
$(error M01 supports CPU MPI only)
endif
endif

# Nombre del ejecutable final
TARGET = grial
ifeq ($(USE_MPI),1)
TARGET = grial_mpi
endif
.DEFAULT_GOAL := all

.PHONY: all grhd2
all: $(TARGET) grhd2

# Compatibilidad con lanzadores anteriores; el ejecutable principal es grial.
grhd2: $(TARGET)
	ln -sfn $(TARGET) $@
TEST_HLLC_TARGET = tests/test_hllc_fluxes
TEST_GW_GEOMETRY_TARGET = tests/test_gw_geometry
TEST_GW_STRESS_TARGET = tests/test_gw_stress
TEST_PPI_GW_TARGET = tests/test_ppi_gw_diagnostics
TEST_FD_TARGET = tests/test_fd_initial
TEST_TORUS_PARAMETERS_TARGET = tests/test_torus_parameters
TEST_FM_TARGET = tests/test_fm_initial
TEST_POLAR_TARGET = tests/test_polar_boundaries
TEST_CFL_TARGET = tests/test_multidimensional_cfl
TEST_MINKOWSKI_SPHERICAL_TARGET = tests/test_minkowski_spherical
TEST_CASE_CORRECTIONS_TARGET = tests/test_case_corrections
TEST_CHECKPOINT_TARGET = tests/test_checkpoint_v5
TEST_OUTPUT_TARGET = tests/test_output_paths
TEST_ACC_RECON_TARGET = tests/test_accelerator_reconstruction
TEST_ACC_RHS_TARGET = tests/test_accelerator_rhs
TEST_NUMERICAL_TARGET = tests/test_numerical_diagnostics
TEST_SPATIAL_TARGET = tests/test_spatial_order
TEST_SENSOR_BOUNDARY_TARGET = tests/test_sensor_boundaries
TEST_SENSOR_SYMMETRY_TARGET = tests/test_sensor_symmetry
TEST_KHI_DIAGNOSTICS_TARGET = tests/test_khi_diagnostics
TEST_MICHEL_INITIAL_TARGET = tests/test_michel_initial
TEST_PRIMITIVE_RECOVERY_TARGET = tests/test_primitive_recovery
TEST_SRHD_CONSERVATION_TARGET = tests/test_srhd_conservation

.PHONY: server-broadwell openacc-p100 openacc-rtx4070 benchmark-jet-p100 test test-hllc test-gw test-gw-geometry test-gw-stress test-fd test-fm test-polar test-cfl test-minkowski-spherical test-cases test-checkpoint test-output test-acc-reconstruction test-acc-rhs clean

# Ejecutable transferible para Intel Xeon E5 v4 (Broadwell-EP).
server-broadwell:
	$(MAKE) -f $(firstword $(MAKEFILE_LIST)) clean
	$(MAKE) -f $(firstword $(MAKEFILE_LIST)) ARCH_FLAGS='-march=broadwell -mtune=broadwell' all

# NVIDIA Tesla P100 = Pascal, compute capability 6.0. No se usa -fast ni
# precisiones mixtas: todos los real*8 permanecen en IEEE FP64.
openacc-p100:
	$(MAKE) -f $(firstword $(MAKEFILE_LIST)) clean
	$(MAKE) -f $(firstword $(MAKEFILE_LIST)) FC=nvfortran FFLAGS='-O3 -acc=gpu -gpu=cc60,lineinfo -Minfo=accel -Mfree -Mextend -cpp -DUSE_OPENACC' all

# Compilación funcional para Zotz/RTX 4070 (Ada, CC 8.9). Conserva FP64;
# sirve para validar kernels, pero no representa el rendimiento FP64 de P100.
openacc-rtx4070:
	$(MAKE) -f $(firstword $(MAKEFILE_LIST)) clean
	$(MAKE) -f $(firstword $(MAKEFILE_LIST)) FC=nvfortran FFLAGS='-O3 -acc=gpu -gpu=cc89,lineinfo -Minfo=accel -Mfree -Mextend -cpp -DUSE_OPENACC' all

# Compila CPU y P100, ejecuta el mismo jet y compara tiempos y campos finales.
benchmark-jet-p100:
	bash benchmarks/jet/run_cpu_p100.sh

# Lista de archivos objeto en el orden estricto de dependencias
EVOLUTION_PARALLEL_OBJS = parallel_runtime.o parallel_timing.o parallel_diagnostics.o
DIAGNOSTIC_OBJS = parallel_runtime.o parallel_timing.o domain_decomposition.o parallel_diagnostics.o
OBJS = parallel_runtime.o parallel_timing.o domain_decomposition.o \
       variables.o parallel_diagnostics.o \
       parameters.o \
       metrics.o \
       output.o \
       initialization.o \
       conditions.o \
       halo_exchange.o \
       equations.o \
       fluxes.o \
       reconstruction.o \
       evolution.o \
       time_integrator.o \
       parallel_perturbations.o \
       parallel_checkpoint.o \
       accelerator.o \
       grial.o

# Regla principal para construir el ejecutable
APPLICATION_OBJS = parallel_observer.o parallel_visualization.o parallel_application.o
$(TARGET): $(OBJS) $(APPLICATION_OBJS)
	$(FC) $(FFLAGS) -o $(TARGET) $(filter-out grial.o,$(OBJS)) $(APPLICATION_OBJS) grial.o

.PHONY: srhd-long-driver test-srhd-long
OBSERVER_OBJS = $(filter-out grial.o,$(OBJS)) parallel_observer.o
VISUALIZATION_OBJS = $(OBSERVER_OBJS) parallel_visualization.o
.PHONY: test-application test-application-mpi
test-application: $(TARGET)
	$(PYTHON) $(SRC_DIR)/tests/test_application_mpi.py --exe ./$(TARGET) --smoke $(if $(filter 0,$(USE_MPI)),--serial)
test: test-application
test-application-mpi: $(TARGET) mpi-segment-driver
ifeq ($(USE_MPI),1)
	$(PYTHON) $(SRC_DIR)/tests/test_application_mpi.py --exe ./$(TARGET) --driver ./tools/mpi_segment_driver $(if $(APPLICATION_OUTPUT),--output $(APPLICATION_OUTPUT))
else
	$(error test-application-mpi requires USE_MPI=1)
endif
.PHONY: mpi-segment-driver
mpi-segment-driver: $(VISUALIZATION_OBJS) $(SRC_DIR)/tools/mpi_segment_driver.f90
	@mkdir -p tools
	$(FC) $(FFLAGS) -o tools/mpi_segment_driver $(SRC_DIR)/tools/mpi_segment_driver.f90 $(VISUALIZATION_OBJS)

VISUALIZATION_TEST = tests/test_parallel_visualization
.PHONY: visualization-driver test-visualization test-visualization-mpi
visualization-driver: $(VISUALIZATION_OBJS) $(SRC_DIR)/tests/test_parallel_visualization.f90
	@mkdir -p tests
	$(FC) $(FFLAGS) -o $(VISUALIZATION_TEST) $(SRC_DIR)/tests/test_parallel_visualization.f90 $(VISUALIZATION_OBJS)
test-visualization: visualization-driver
	@fixture_dir=$$(mktemp -d /tmp/grial_visualization_XXXXXX); cd "$$fixture_dir" && OMP_NUM_THREADS=2 $(abspath $(VISUALIZATION_TEST)) 1 1 4 normal && $(PYTHON) $(abspath $(SRC_DIR))/tools/parallel_vtk.py "$$fixture_dir/good" > "$$fixture_dir/inspection.json"
test: test-visualization
test-visualization-mpi: visualization-driver mpi-segment-driver
ifeq ($(USE_MPI),1)
	$(PYTHON) $(SRC_DIR)/tests/test_visualization_mpi.py --probe ./$(VISUALIZATION_TEST) --driver ./tools/mpi_segment_driver $(if $(VISUALIZATION_OUTPUT),--output $(VISUALIZATION_OUTPUT)) $(VISUALIZATION_FLAGS)
else
	$(error test-visualization-mpi requires USE_MPI=1)
endif

OBSERVER_TEST = tests/test_parallel_observer
.PHONY: observer-driver test-observer test-observer-mpi
observer-driver: $(OBSERVER_OBJS) $(SRC_DIR)/tests/test_parallel_observer.f90
	@mkdir -p tests
	$(FC) $(FFLAGS) -o $(OBSERVER_TEST) $(SRC_DIR)/tests/test_parallel_observer.f90 $(OBSERVER_OBJS)
test-observer: observer-driver
	@fixture_dir=$$(mktemp -d /tmp/grial_observer_XXXXXX); cd "$$fixture_dir" && OMP_NUM_THREADS=2 $(abspath $(OBSERVER_TEST)) 1 1 roundtrip
test: test-observer
.PHONY: test-mpi-segments
test-mpi-segments:
	$(PYTHON) $(SRC_DIR)/tests/test_mpi_segments.py
test: test-mpi-segments
test-observer-mpi: mpi-segment-driver observer-driver
ifeq ($(USE_MPI),1)
	$(PYTHON) $(SRC_DIR)/tests/test_observer_mpi.py --driver ./tools/mpi_segment_driver --probe ./$(OBSERVER_TEST) $(if $(OBSERVER_OUTPUT),--output $(OBSERVER_OUTPUT)) $(OBSERVER_FLAGS)
else
	$(error test-observer-mpi requires USE_MPI=1)
endif

test: test-checkpoint-parallel
srhd-long-driver: $(filter-out grial.o,$(OBJS)) $(SRC_DIR)/tools/srhd_long_driver.f90
	@mkdir -p tools
	$(FC) $(FFLAGS) -o tools/srhd_long_driver $(SRC_DIR)/tools/srhd_long_driver.f90 $(filter-out grial.o,$(OBJS))

test-srhd-long:
	$(PYTHON) $(SRC_DIR)/tests/test_srhd_long.py

test: test-perturbations test-srhd-long test-observables test-diagnostics test-rtheta test-rk test-time test-interfaces test-halos test-domain test-parallel test-segments test-walltime test-hllc test-gw test-fd test-fm test-torus-parameters test-polar test-cfl test-minkowski-spherical test-cases test-checkpoint test-output test-acc-reconstruction test-acc-rhs test-numerical-diagnostics test-spatial-order test-sensor-boundaries test-sensor-symmetry test-khi-diagnostics test-michel-initial test-primitive-recovery test-srhd-conservation

.PHONY: test-srhd-conservation
test-srhd-conservation: variables.o parameters.o metrics.o conditions.o equations.o fluxes.o reconstruction.o evolution.o $(EVOLUTION_PARALLEL_OBJS) initialization.o domain_decomposition.o $(SRC_DIR)/tests/test_srhd_conservation.f90
	$(FC) $(FFLAGS) -o $(TEST_SRHD_CONSERVATION_TARGET) $(SRC_DIR)/tests/test_srhd_conservation.f90 variables.o parameters.o metrics.o conditions.o equations.o fluxes.o reconstruction.o evolution.o $(EVOLUTION_PARALLEL_OBJS) initialization.o domain_decomposition.o
	OMP_NUM_THREADS=1 ./$(TEST_SRHD_CONSERVATION_TARGET)

.PHONY: test-primitive-recovery
test-primitive-recovery: variables.o metrics.o equations.o $(SRC_DIR)/tests/test_primitive_recovery.f90
	$(FC) $(FFLAGS) -o $(TEST_PRIMITIVE_RECOVERY_TARGET) $(SRC_DIR)/tests/test_primitive_recovery.f90 variables.o metrics.o equations.o
	OMP_NUM_THREADS=1 ./$(TEST_PRIMITIVE_RECOVERY_TARGET)

.PHONY: test-sensor-symmetry
test-sensor-symmetry: variables.o reconstruction.o $(SRC_DIR)/tests/test_sensor_symmetry.f90
	$(FC) $(FFLAGS) -o $(TEST_SENSOR_SYMMETRY_TARGET) $(SRC_DIR)/tests/test_sensor_symmetry.f90 variables.o reconstruction.o
	OMP_NUM_THREADS=4 ./$(TEST_SENSOR_SYMMETRY_TARGET)

.PHONY: test-michel-initial
test-michel-initial: variables.o metrics.o conditions.o $(SRC_DIR)/tests/test_michel_initial.f90
	$(FC) $(FFLAGS) -o $(TEST_MICHEL_INITIAL_TARGET) $(SRC_DIR)/tests/test_michel_initial.f90 variables.o metrics.o conditions.o
	OMP_NUM_THREADS=1 ./$(TEST_MICHEL_INITIAL_TARGET)

.PHONY: test-khi-diagnostics
test-khi-diagnostics: variables.o metrics.o conditions.o output.o $(DIAGNOSTIC_OBJS) $(SRC_DIR)/tests/test_khi_diagnostics.f90
	$(FC) $(FFLAGS) -o $(TEST_KHI_DIAGNOSTICS_TARGET) $(SRC_DIR)/tests/test_khi_diagnostics.f90 variables.o metrics.o conditions.o output.o $(DIAGNOSTIC_OBJS)
	OMP_NUM_THREADS=1 ./$(TEST_KHI_DIAGNOSTICS_TARGET)

.PHONY: test-sensor-boundaries
test-sensor-boundaries: variables.o parameters.o parallel_runtime.o metrics.o conditions.o equations.o fluxes.o reconstruction.o initialization.o domain_decomposition.o $(SRC_DIR)/tests/test_sensor_boundaries.f90
	$(FC) $(FFLAGS) -o $(TEST_SENSOR_BOUNDARY_TARGET) $(SRC_DIR)/tests/test_sensor_boundaries.f90 variables.o parameters.o parallel_runtime.o metrics.o conditions.o equations.o fluxes.o reconstruction.o initialization.o domain_decomposition.o
	OMP_NUM_THREADS=1 ./$(TEST_SENSOR_BOUNDARY_TARGET)

.PHONY: test-spatial-order
test-spatial-order: variables.o metrics.o equations.o fluxes.o reconstruction.o $(SRC_DIR)/tests/test_spatial_order.f90
	$(FC) $(FFLAGS) -o $(TEST_SPATIAL_TARGET) $(SRC_DIR)/tests/test_spatial_order.f90 variables.o metrics.o equations.o fluxes.o reconstruction.o
	OMP_NUM_THREADS=1 ./$(TEST_SPATIAL_TARGET)

.PHONY: test-numerical-diagnostics
test-numerical-diagnostics: variables.o metrics.o equations.o fluxes.o reconstruction.o $(SRC_DIR)/tests/test_numerical_diagnostics.f90
	$(FC) $(FFLAGS) -o $(TEST_NUMERICAL_TARGET) $(SRC_DIR)/tests/test_numerical_diagnostics.f90 variables.o metrics.o equations.o fluxes.o reconstruction.o
	OMP_NUM_THREADS=4 ./$(TEST_NUMERICAL_TARGET)

test-hllc: variables.o metrics.o equations.o fluxes.o $(SRC_DIR)/tests/test_hllc_fluxes.f90
	$(FC) $(FFLAGS) -o $(TEST_HLLC_TARGET) $(SRC_DIR)/tests/test_hllc_fluxes.f90 variables.o metrics.o equations.o fluxes.o
	./$(TEST_HLLC_TARGET)

test-gw-geometry: variables.o metrics.o conditions.o output.o $(DIAGNOSTIC_OBJS) $(SRC_DIR)/tests/test_gw_geometry.f90
	$(FC) $(FFLAGS) -o $(TEST_GW_GEOMETRY_TARGET) $(SRC_DIR)/tests/test_gw_geometry.f90 variables.o metrics.o conditions.o output.o $(DIAGNOSTIC_OBJS)
	./$(TEST_GW_GEOMETRY_TARGET)

test-gw-stress: variables.o metrics.o conditions.o output.o $(DIAGNOSTIC_OBJS) $(SRC_DIR)/tests/test_gw_stress.f90
	$(FC) $(FFLAGS) -o $(TEST_GW_STRESS_TARGET) $(SRC_DIR)/tests/test_gw_stress.f90 variables.o metrics.o conditions.o output.o $(DIAGNOSTIC_OBJS)
	./$(TEST_GW_STRESS_TARGET)

test-gw: test-gw-geometry test-gw-stress test-ppi-gw

.PHONY: test-ppi-gw
test-ppi-gw: variables.o metrics.o conditions.o equations.o output.o $(DIAGNOSTIC_OBJS) $(SRC_DIR)/tests/test_ppi_gw_diagnostics.f90
	$(FC) $(FFLAGS) -o $(TEST_PPI_GW_TARGET) $(SRC_DIR)/tests/test_ppi_gw_diagnostics.f90 variables.o metrics.o conditions.o equations.o output.o $(DIAGNOSTIC_OBJS)
	@fixture_dir=$$(mktemp -d /tmp/grial_ppi_gw_XXXXXX); OMP_NUM_THREADS=1 ./$(TEST_PPI_GW_TARGET) "$$fixture_dir"

test-fd: variables.o metrics.o conditions.o $(SRC_DIR)/tests/test_fd_initial.f90
	$(FC) $(FFLAGS) -o $(TEST_FD_TARGET) $(SRC_DIR)/tests/test_fd_initial.f90 variables.o metrics.o conditions.o
	./$(TEST_FD_TARGET)

.PHONY: torus-profiles
torus-profiles: variables.o metrics.o conditions.o $(SRC_DIR)/analysis/torus_initial_profiles.f90
	$(FC) $(FFLAGS) -o analysis/torus_initial_profiles $(SRC_DIR)/analysis/torus_initial_profiles.f90 variables.o metrics.o conditions.o

.PHONY: test-torus-parameters
test-torus-parameters: variables.o parameters.o parallel_runtime.o $(SRC_DIR)/tests/test_torus_parameters.f90
	$(FC) $(FFLAGS) -o $(TEST_TORUS_PARAMETERS_TARGET) $(SRC_DIR)/tests/test_torus_parameters.f90 variables.o parameters.o parallel_runtime.o
	./$(TEST_TORUS_PARAMETERS_TARGET)

test-fm: variables.o metrics.o conditions.o $(SRC_DIR)/tests/test_fm_initial.f90
	$(FC) $(FFLAGS) -o $(TEST_FM_TARGET) $(SRC_DIR)/tests/test_fm_initial.f90 variables.o metrics.o conditions.o
	./$(TEST_FM_TARGET)

test-polar: variables.o metrics.o conditions.o $(SRC_DIR)/tests/test_polar_boundaries.f90
	$(FC) $(FFLAGS) -o $(TEST_POLAR_TARGET) $(SRC_DIR)/tests/test_polar_boundaries.f90 variables.o metrics.o conditions.o
	./$(TEST_POLAR_TARGET)

test-cfl: variables.o metrics.o equations.o fluxes.o reconstruction.o evolution.o $(EVOLUTION_PARALLEL_OBJS) domain_decomposition.o $(SRC_DIR)/tests/test_multidimensional_cfl.f90
	$(FC) $(FFLAGS) -o $(TEST_CFL_TARGET) $(SRC_DIR)/tests/test_multidimensional_cfl.f90 variables.o metrics.o equations.o fluxes.o reconstruction.o evolution.o $(EVOLUTION_PARALLEL_OBJS) domain_decomposition.o
	./$(TEST_CFL_TARGET)

test-minkowski-spherical: variables.o metrics.o $(SRC_DIR)/tests/test_minkowski_spherical.f90
	$(FC) $(FFLAGS) -o $(TEST_MINKOWSKI_SPHERICAL_TARGET) $(SRC_DIR)/tests/test_minkowski_spherical.f90 variables.o metrics.o
	./$(TEST_MINKOWSKI_SPHERICAL_TARGET)

test-cases: variables.o parameters.o parallel_runtime.o metrics.o conditions.o equations.o fluxes.o initialization.o domain_decomposition.o $(SRC_DIR)/tests/test_case_corrections.f90
	$(FC) $(FFLAGS) -o $(TEST_CASE_CORRECTIONS_TARGET) $(SRC_DIR)/tests/test_case_corrections.f90 variables.o parameters.o parallel_runtime.o metrics.o conditions.o equations.o fluxes.o initialization.o domain_decomposition.o
	./$(TEST_CASE_CORRECTIONS_TARGET)

test-checkpoint: variables.o metrics.o conditions.o output.o $(DIAGNOSTIC_OBJS) $(SRC_DIR)/tests/test_checkpoint_v5.f90
	$(FC) $(FFLAGS) -o $(TEST_CHECKPOINT_TARGET) $(SRC_DIR)/tests/test_checkpoint_v5.f90 variables.o metrics.o conditions.o output.o $(DIAGNOSTIC_OBJS)
	./$(TEST_CHECKPOINT_TARGET)

test-output: variables.o metrics.o conditions.o output.o $(DIAGNOSTIC_OBJS) $(SRC_DIR)/tests/test_output_paths.f90
	$(FC) $(FFLAGS) -o $(TEST_OUTPUT_TARGET) $(SRC_DIR)/tests/test_output_paths.f90 variables.o metrics.o conditions.o output.o $(DIAGNOSTIC_OBJS)
	./$(TEST_OUTPUT_TARGET)

test-acc-reconstruction: variables.o reconstruction.o $(SRC_DIR)/tests/test_accelerator_reconstruction.f90
	$(FC) $(FFLAGS) -o $(TEST_ACC_RECON_TARGET) $(SRC_DIR)/tests/test_accelerator_reconstruction.f90 variables.o reconstruction.o
	./$(TEST_ACC_RECON_TARGET)

test-acc-rhs: variables.o metrics.o equations.o fluxes.o reconstruction.o evolution.o $(EVOLUTION_PARALLEL_OBJS) domain_decomposition.o $(SRC_DIR)/tests/test_accelerator_rhs.f90
	$(FC) $(FFLAGS) -o $(TEST_ACC_RHS_TARGET) $(SRC_DIR)/tests/test_accelerator_rhs.f90 variables.o metrics.o equations.o fluxes.o reconstruction.o evolution.o $(EVOLUTION_PARALLEL_OBJS) domain_decomposition.o
	./$(TEST_ACC_RHS_TARGET)

# Regla genérica para compilar archivos .f90 a objetos .o
%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c $<

# --- DEPENDENCIAS DE LOS MÓDULOS ---
# Esto le dice a Make qué debe compilarse antes que otra cosa
parallel_timing.o: parallel_runtime.o
metrics.o: variables.o
parameters.o: variables.o parallel_runtime.o
output.o: variables.o metrics.o conditions.o $(DIAGNOSTIC_OBJS)
initialization.o: parallel_runtime.o domain_decomposition.o variables.o parameters.o metrics.o conditions.o fluxes.o
halo_exchange.o: parallel_runtime.o domain_decomposition.o variables.o conditions.o
conditions.o: variables.o metrics.o
equations.o: variables.o metrics.o
fluxes.o: variables.o equations.o metrics.o
reconstruction.o: variables.o
parallel_diagnostics.o: parallel_runtime.o parallel_timing.o domain_decomposition.o variables.o
evolution.o: parallel_diagnostics.o parallel_timing.o variables.o reconstruction.o fluxes.o equations.o
accelerator.o: variables.o
time_integrator.o: parallel_diagnostics.o variables.o equations.o evolution.o conditions.o output.o halo_exchange.o parallel_timing.o
parallel_checkpoint.o: parallel_runtime.o variables.o domain_decomposition.o parallel_diagnostics.o halo_exchange.o output.o parameters.o initialization.o metrics.o equations.o
parallel_observer.o: parallel_checkpoint.o parallel_diagnostics.o output.o
parallel_visualization.o: parallel_runtime.o domain_decomposition.o variables.o parallel_diagnostics.o output.o
parallel_perturbations.o: parallel_runtime.o parallel_timing.o parallel_diagnostics.o domain_decomposition.o variables.o conditions.o equations.o halo_exchange.o output.o
parallel_application.o: parallel_observer.o parallel_visualization.o parallel_perturbations.o initialization.o halo_exchange.o time_integrator.o accelerator.o

grial.o: parallel_application.o time_integrator.o parallel_timing.o parallel_runtime.o variables.o parameters.o metrics.o initialization.o conditions.o equations.o fluxes.o reconstruction.o evolution.o accelerator.o output.o

# Limpiar archivos compilados (comando: make clean)
clean:
	rm -f $(VISUALIZATION_TEST)
	rm -f $(OBSERVER_TEST) tools/mpi_segment_driver
	rm -f $(PERTURBATIONS_TEST) $(PARALLEL_CHECKPOINT_TEST) $(PARALLEL_RESTART_TEST)
	rm -f tools/srhd_long_driver
	rm -f $(OBSERVABLES_TEST) $(DIAGNOSTICS_TEST) $(RTHETA_TEST) $(RK_TEST) $(TIMING_TEST) $(INTERFACE_TEST) $(HALO_TEST) $(DOMAIN_TEST) $(PARALLEL_TEST)
	rm -f $(TEST_PPI_GW_TARGET)
	rm -f analysis/torus_initial_profiles analysis/mpi_shock_profiles
	rm -f *.o *.mod $(TARGET) grhd2 $(TEST_HLLC_TARGET) $(TEST_GW_GEOMETRY_TARGET) $(TEST_GW_STRESS_TARGET) $(TEST_FD_TARGET) $(TEST_FM_TARGET) $(TEST_TORUS_PARAMETERS_TARGET) $(TEST_POLAR_TARGET) $(TEST_CFL_TARGET) $(TEST_MINKOWSKI_SPHERICAL_TARGET) $(TEST_CASE_CORRECTIONS_TARGET) $(TEST_CHECKPOINT_TARGET) $(TEST_OUTPUT_TARGET) $(TEST_ACC_RECON_TARGET) $(TEST_ACC_RHS_TARGET) $(TEST_NUMERICAL_TARGET) $(TEST_SPATIAL_TARGET) $(TEST_SENSOR_BOUNDARY_TARGET) $(TEST_SENSOR_SYMMETRY_TARGET) $(TEST_KHI_DIAGNOSTICS_TARGET) $(TEST_MICHEL_INITIAL_TARGET) $(TEST_PRIMITIVE_RECOVERY_TARGET) $(TEST_SRHD_CONSERVATION_TARGET)

.PHONY: test-walltime
test-walltime: $(TARGET)
	$(PYTHON) $(SRC_DIR)/tests/test_walltime_restart.py --exe ./$(TARGET) $(if $(filter 1,$(USE_MPI)),--mpi-one)

.PHONY: test-segments
test-segments: $(TARGET) test-checkpoint
	$(PYTHON) $(SRC_DIR)/tests/test_segment_launcher.py --exe ./$(TARGET) --probe ./$(TEST_CHECKPOINT_TARGET)

.PHONY: test-parallel
PARALLEL_TEST = tests/test_parallel_runtime
test-parallel: parallel_runtime.o $(SRC_DIR)/tests/test_parallel_runtime.f90
	$(FC) $(FFLAGS) -o $(PARALLEL_TEST) $(SRC_DIR)/tests/test_parallel_runtime.f90 parallel_runtime.o
	OMP_NUM_THREADS=2 ./$(PARALLEL_TEST) ok

.PHONY: test-domain test-domain-mpi
DOMAIN_TEST = tests/test_domain_decomposition
test-domain: parallel_runtime.o domain_decomposition.o variables.o parameters.o metrics.o conditions.o equations.o fluxes.o initialization.o $(SRC_DIR)/tests/test_domain_decomposition.f90
	$(FC) $(FFLAGS) -o $(DOMAIN_TEST) $(SRC_DIR)/tests/test_domain_decomposition.f90 parallel_runtime.o domain_decomposition.o variables.o parameters.o metrics.o conditions.o equations.o fluxes.o initialization.o
	OMP_NUM_THREADS=2 ./$(DOMAIN_TEST)

test-domain-mpi: test-domain
ifeq ($(USE_MPI),1)
	$(PYTHON) $(SRC_DIR)/tests/test_domain_mpi.py --probe ./$(DOMAIN_TEST)
else
	$(error test-domain-mpi requires USE_MPI=1)
endif

.PHONY: test-halos test-halos-mpi
HALO_TEST = tests/test_halo_exchange
test-halos: parallel_runtime.o domain_decomposition.o halo_exchange.o variables.o parameters.o metrics.o conditions.o equations.o fluxes.o initialization.o $(SRC_DIR)/tests/test_halo_exchange.f90
	$(FC) $(FFLAGS) -o $(HALO_TEST) $(SRC_DIR)/tests/test_halo_exchange.f90 parallel_runtime.o domain_decomposition.o halo_exchange.o variables.o parameters.o metrics.o conditions.o equations.o fluxes.o initialization.o
	OMP_NUM_THREADS=2 ./$(HALO_TEST)

test-halos-mpi: test-halos
ifeq ($(USE_MPI),1)
	$(PYTHON) $(SRC_DIR)/tests/test_halo_mpi.py --probe ./$(HALO_TEST)
else
	$(error test-halos-mpi requires USE_MPI=1)
endif

.PHONY: test-interfaces test-interfaces-mpi
INTERFACE_TEST = tests/test_mpi_interfaces
test-interfaces: $(filter-out grial.o time_integrator.o parallel_perturbations.o parallel_checkpoint.o accelerator.o output.o,$(OBJS)) $(SRC_DIR)/tests/test_mpi_interfaces.f90
	$(FC) $(FFLAGS) -o $(INTERFACE_TEST) $(SRC_DIR)/tests/test_mpi_interfaces.f90 $(filter-out grial.o time_integrator.o parallel_perturbations.o parallel_checkpoint.o accelerator.o output.o,$(OBJS))
	OMP_NUM_THREADS=2 ./$(INTERFACE_TEST)

test-interfaces-mpi: test-interfaces
ifeq ($(USE_MPI),1)
	$(PYTHON) $(SRC_DIR)/tests/test_interfaces_mpi.py --probe ./$(INTERFACE_TEST)
else
	$(error test-interfaces-mpi requires USE_MPI=1)
endif

.PHONY: test-observables test-observables-mpi
PARALLEL_RESTART_TEST = tests/test_parallel_restart
.PHONY: restart-driver test-restart-parallel test-restart-mpi
test: test-restart-parallel
restart-driver: $(filter-out grial.o accelerator.o,$(OBJS)) $(SRC_DIR)/tests/test_parallel_restart.f90
	$(FC) $(FFLAGS) -o $(PARALLEL_RESTART_TEST) $(SRC_DIR)/tests/test_parallel_restart.f90 $(filter-out grial.o accelerator.o,$(OBJS))
test-restart-parallel: restart-driver
	$(PYTHON) $(SRC_DIR)/tests/test_restart_mpi.py --probe $(abspath $(PARALLEL_RESTART_TEST)) --serial --smoke

test-restart-mpi: test-restart-parallel
ifeq ($(USE_MPI),1)
	$(PYTHON) $(SRC_DIR)/tests/test_restart_mpi.py --probe $(abspath $(PARALLEL_RESTART_TEST)) $(if $(RESTART_MPI_OUTPUT),--output $(RESTART_MPI_OUTPUT))
else
	$(error test-restart-mpi requires USE_MPI=1)
endif

PARALLEL_CHECKPOINT_TEST = tests/test_parallel_checkpoint
.PHONY: test-checkpoint-parallel test-checkpoint-mpi
test-checkpoint-parallel: $(filter-out grial.o accelerator.o,$(OBJS)) $(SRC_DIR)/tests/test_parallel_checkpoint.f90
	$(FC) $(FFLAGS) -o $(PARALLEL_CHECKPOINT_TEST) $(SRC_DIR)/tests/test_parallel_checkpoint.f90 $(filter-out grial.o accelerator.o,$(OBJS))
	@fixture_dir=$$(mktemp -d /tmp/grial_parallel_checkpoint_XXXXXX); cd "$$fixture_dir" && OMP_NUM_THREADS=2 $(abspath $(PARALLEL_CHECKPOINT_TEST)) 1 1 report

test-checkpoint-mpi: test-checkpoint-parallel
ifeq ($(USE_MPI),1)
	$(PYTHON) $(SRC_DIR)/tests/test_checkpoint_mpi.py --probe ./$(PARALLEL_CHECKPOINT_TEST) $(if $(CHECKPOINT_MPI_OUTPUT),--output $(CHECKPOINT_MPI_OUTPUT))
else
	$(error test-checkpoint-mpi requires USE_MPI=1)
endif

PERTURBATIONS_TEST = tests/test_parallel_perturbations
.PHONY: test-perturbations test-perturbations-mpi
test-perturbations: $(filter-out grial.o accelerator.o,$(OBJS)) $(SRC_DIR)/tests/test_parallel_perturbations.f90
	$(FC) $(FFLAGS) -o $(PERTURBATIONS_TEST) $(SRC_DIR)/tests/test_parallel_perturbations.f90 $(filter-out grial.o accelerator.o,$(OBJS))
	@fixture_dir=$$(mktemp -d /tmp/grial_perturbations_XXXXXX); cd "$$fixture_dir" && OMP_NUM_THREADS=2 $(abspath $(PERTURBATIONS_TEST)) 1 1 report

test-perturbations-mpi: test-perturbations
ifeq ($(USE_MPI),1)
	$(PYTHON) $(SRC_DIR)/tests/test_perturbations_mpi.py --probe ./$(PERTURBATIONS_TEST) $(if $(PERTURBATIONS_OUTPUT),--output $(PERTURBATIONS_OUTPUT))
else
	$(error test-perturbations-mpi requires USE_MPI=1)
endif

OBSERVABLES_TEST = tests/test_parallel_observables
test-observables: $(filter-out grial.o accelerator.o,$(OBJS)) $(SRC_DIR)/tests/test_parallel_observables.f90
	$(FC) $(FFLAGS) -o $(OBSERVABLES_TEST) $(SRC_DIR)/tests/test_parallel_observables.f90 $(filter-out grial.o accelerator.o,$(OBJS))
	@fixture_dir=$$(mktemp -d /tmp/grial_observables_XXXXXX); cd "$$fixture_dir" && OMP_NUM_THREADS=2 $(abspath $(OBSERVABLES_TEST)) 1 1 report

test-observables-mpi: test-observables
ifeq ($(USE_MPI),1)
	$(PYTHON) $(SRC_DIR)/tests/test_observables_mpi.py --probe ./$(OBSERVABLES_TEST) $(if $(OBSERVABLES_OUTPUT),--output $(OBSERVABLES_OUTPUT)) $(OBSERVABLES_FLAGS)
else
	$(error test-observables-mpi requires USE_MPI=1)
endif

.PHONY: test-diagnostics test-diagnostics-mpi
DIAGNOSTICS_TEST = tests/test_parallel_diagnostics
test-diagnostics: $(filter-out grial.o accelerator.o,$(OBJS)) $(SRC_DIR)/tests/test_parallel_diagnostics.f90
	$(FC) $(FFLAGS) -o $(DIAGNOSTICS_TEST) $(SRC_DIR)/tests/test_parallel_diagnostics.f90 $(filter-out grial.o accelerator.o,$(OBJS))
	OMP_NUM_THREADS=2 ./$(DIAGNOSTICS_TEST)

test-diagnostics-mpi: test-diagnostics
ifeq ($(USE_MPI),1)
	$(PYTHON) $(SRC_DIR)/tests/test_diagnostics_mpi.py --probe ./$(DIAGNOSTICS_TEST) $(if $(DIAGNOSTICS_OUTPUT),--output $(DIAGNOSTICS_OUTPUT)) $(DIAGNOSTICS_FLAGS)
else
	$(error test-diagnostics-mpi requires USE_MPI=1)
endif

.PHONY: test-rtheta test-rtheta-mpi
RTHETA_TEST = tests/test_rtheta
test-rtheta: $(filter-out grial.o accelerator.o,$(OBJS)) $(SRC_DIR)/tests/test_rtheta.f90
	$(FC) $(FFLAGS) -o $(RTHETA_TEST) $(SRC_DIR)/tests/test_rtheta.f90 $(filter-out grial.o accelerator.o,$(OBJS))
	OMP_NUM_THREADS=2 ./$(RTHETA_TEST)

test-rtheta-mpi: test-rtheta
ifeq ($(USE_MPI),1)
	$(PYTHON) $(SRC_DIR)/tests/test_rtheta_mpi.py --probe ./$(RTHETA_TEST) $(if $(RTHETA_OUTPUT),--output $(RTHETA_OUTPUT)) $(RTHETA_FLAGS)
else
	$(error test-rtheta-mpi requires USE_MPI=1)
endif

.PHONY: test-rk test-rk-mpi
analysis/mpi_shock_profiles: $(filter-out grial.o accelerator.o,$(OBJS)) $(SRC_DIR)/analysis/mpi_shock_profiles.f90
	$(FC) $(FFLAGS) -o $@ $(SRC_DIR)/analysis/mpi_shock_profiles.f90 $(filter-out grial.o accelerator.o,$(OBJS))

RK_TEST = tests/test_distributed_rk
test-rk: $(filter-out grial.o accelerator.o,$(OBJS)) $(SRC_DIR)/tests/test_distributed_rk.f90
	$(FC) $(FFLAGS) -o $(RK_TEST) $(SRC_DIR)/tests/test_distributed_rk.f90 $(filter-out grial.o accelerator.o,$(OBJS))
	OMP_NUM_THREADS=2 ./$(RK_TEST)

test-rk-mpi: test-rk
ifeq ($(USE_MPI),1)
	$(PYTHON) $(SRC_DIR)/tests/test_rk_mpi.py --probe ./$(RK_TEST)
else
	$(error test-rk-mpi requires USE_MPI=1)
endif

.PHONY: test-time test-time-mpi
TIMING_TEST = tests/test_parallel_timing
test-time: $(filter-out grial.o time_integrator.o parallel_perturbations.o parallel_checkpoint.o accelerator.o output.o,$(OBJS)) $(SRC_DIR)/tests/test_parallel_timing.f90
	$(FC) $(FFLAGS) -o $(TIMING_TEST) $(SRC_DIR)/tests/test_parallel_timing.f90 $(filter-out grial.o time_integrator.o parallel_perturbations.o parallel_checkpoint.o accelerator.o output.o,$(OBJS))
	OMP_NUM_THREADS=2 ./$(TIMING_TEST)

test-time-mpi: test-time
ifeq ($(USE_MPI),1)
	$(PYTHON) $(SRC_DIR)/tests/test_timing_mpi.py --probe ./$(TIMING_TEST)
else
	$(error test-time-mpi requires USE_MPI=1)
endif

.PHONY: test-mpi
test-mpi: $(TARGET) test-parallel
ifeq ($(USE_MPI),1)
	$(PYTHON) $(SRC_DIR)/tests/test_mpi_optional.py --exe ./$(TARGET) --probe ./$(PARALLEL_TEST) $(if $(wildcard $(SRC_DIR)/grial),--serial-exe $(SRC_DIR)/grial,)
else
	$(error test-mpi requires USE_MPI=1)
endif

endif
