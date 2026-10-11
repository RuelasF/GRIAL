#!/usr/bin/env bash
#SBATCH --job-name=ppi-w3-smoke
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --array=0-2%3
#SBATCH --output=slurm_logs/weno3-smoke-%A_%a.log
#SBATCH --mail-user=francisco.ruelas9565@alumnos.udg.mx
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS

set -euo pipefail
umask 027

# Igual que los trabajos probados de OllinSphere: se envia desde la raiz del
# repositorio y se usa un ejecutable compilado previamente en la sesion ollin.
REPO=${SLURM_SUBMIT_DIR:?Envie este sbatch desde la raiz del repositorio GRHD}
test -f "$REPO/Makefile"
test -x "$REPO/grial"

cases=(ef_a00_weno3 ks_a00_weno3 ks_a02_weno3)
label=${cases[$SLURM_ARRAY_TASK_ID]}
case "$label" in
  ef_a00_weno3)
    metric=eddington_finkelstein
    geometry=spherical
    spin=0.0
    ;;
  ks_a00_weno3)
    metric=kerr_schild
    geometry=spheroidal
    spin=0.0
    ;;
  ks_a02_weno3)
    metric=kerr_schild
    geometry=spheroidal
    spin=0.2
    ;;
esac

RUN_ROOT="$REPO/ppi_weno3_smoke/job_${SLURM_ARRAY_JOB_ID}"
CASE_ROOT="$RUN_ROOT/$label"
PARAMETER_ROOT="$RUN_ROOT/parameters"
PARAMETER_FILE="$PARAMETER_ROOT/$label.par"
OUTPUT_ROOT="$CASE_ROOT/output"
PROGRAM_LOG="$CASE_ROOT/grial.log"
mkdir -p "$REPO/slurm_logs" "$PARAMETER_ROOT" "$CASE_ROOT"

cat >"$PARAMETER_FILE" <<EOF
problem = fishbone_equatorial
restart = false
metric = $metric
geometry = $geometry
bh_mass = 1.0
spin = $spin
adiabatic_index = 1.333333333333333333
fm_inner_radius = 6.25
fm_pressure_max_radius = 9.25
fm_polytropic_constant = 0.0015
nx = 80
ny = 1
nz = 40
r_min = 1.2
r_max = 40.0
theta_min = 0.0
theta_max = 1.0
phi_min = 0.0
phi_max = 2.0
logarithmic_r = true
reconstruction = weno3
riemann_solver = hlle
tvd_limiter = mc
shock_sensor = false
final_time = 0.5
cfl = 0.4
save_interval = 0.5
checkpoint_interval = 10.0
output_prefix = $label
output_folder = $OUTPUT_ROOT
vtk_mapping = physical
extract_gw = true
extract_mdot = true
ppi_diagnostics = true
diagnostic_stride = 10
extraction_stride = 10
mass_monitor_stride = 0
perturbation = pressure_noise
perturbation_seed = 3435
perturbation_time = 1000.0
perturbation_amplitude = 0.01
perturbation_mode = 1.0
EOF

"$REPO/grial" --check "$PARAMETER_FILE"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_DYNAMIC=FALSE
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_DISPLAY_ENV=TRUE
export OMP_DISPLAY_AFFINITY=TRUE
export GFORTRAN_UNBUFFERED_ALL=y

printf 'START %s case=%s host=%s task=%s threads=%s cpus=%s\n' \
  "$(date --iso-8601=seconds)" "$label" "$(hostname)" \
  "$SLURM_ARRAY_TASK_ID" "$OMP_NUM_THREADS" \
  "$(awk '/^Cpus_allowed_list:/ {print $2}' /proc/self/status)"

srun --ntasks=1 --cpus-per-task="$OMP_NUM_THREADS" --cpu-bind=cores \
  "$REPO/grial" "$PARAMETER_FILE" >"$PROGRAM_LOG" 2>&1

grep -q 'FINALIZADA EXITOSAMENTE' "$PROGRAM_LOG"
test -s "$OUTPUT_ROOT/weno3_hlle/${label}_weno3_hlle_0001.vtk"
test -s "$OUTPUT_ROOT/weno3_hlle/global_diagnostics.dat"
last_time=$(awk '!/^#/ && NF {value=$1} END {print value}' \
  "$OUTPUT_ROOT/weno3_hlle/global_diagnostics.dat")
awk -v value="$last_time" 'BEGIN {exit(value >= 0.499999999 && value <= 0.500000001 ? 0 : 1)}'

printf 'PASS %s case=%s t=%s output=%s\n' \
  "$(date --iso-8601=seconds)" "$label" "$last_time" "$CASE_ROOT"
