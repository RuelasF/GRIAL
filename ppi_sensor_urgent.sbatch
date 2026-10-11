#!/usr/bin/env bash
#SBATCH --job-name=ppi-sensor
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=16:00:00
#SBATCH --array=0-7%8
#SBATCH --output=slurm_logs/ppi-sensor-%A_%a.log
#SBATCH --mail-user=francisco.ruelas9565@alumnos.udg.mx
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS

set -Eeuo pipefail
umask 027

readonly REPO="${GRHD_REPO_ROOT:-${SLURM_SUBMIT_DIR:?Envíe desde la raíz de GRHD}}"
readonly EXECUTABLE="${GRHD_EXECUTABLE:-$REPO/grial}"
readonly CAMPAIGN_ROOT="${PPI_SENSOR_ROOT:-$REPO/ppi_sensor_mp5_weno5_300x1x150_t2000_cluster}"
readonly TASK_ID="${SLURM_ARRAY_TASK_ID:?Falta SLURM_ARRAY_TASK_ID}"

[[ -x "$EXECUTABLE" ]] || {
  printf 'ERROR: no existe el ejecutable %s\n' "$EXECUTABLE" >&2
  exit 1
}

cases=(
  ef_a00_weno5 ef_a00_mp5
  ks_a00_weno5 ks_a00_mp5
  ks_a05_weno5 ks_a05_mp5
  ks_a09_weno5 ks_a09_mp5
)
(( TASK_ID >= 0 && TASK_ID < ${#cases[@]} )) || {
  printf 'ERROR: índice fuera de rango: %s\n' "$TASK_ID" >&2
  exit 2
}
label=${cases[$TASK_ID]}
reconstruction=${label##*_}

if [[ "$label" == ef_* ]]; then
  metric=eddington_finkelstein
  geometry=spherical
  spin=0.0
else
  metric=kerr_schild
  geometry=spheroidal
  case "$label" in
    ks_a00_*) spin=0.0 ;;
    ks_a05_*) spin=0.5 ;;
    ks_a09_*) spin=0.9 ;;
    *) printf 'ERROR: espín desconocido en %s\n' "$label" >&2; exit 2 ;;
  esac
fi

readonly PARAMETER_ROOT="$CAMPAIGN_ROOT/parameters"
readonly CASE_ROOT="$CAMPAIGN_ROOT/results/$label/session_000"
readonly STATE_ROOT="$CAMPAIGN_ROOT/state"
readonly PARAMETER_FILE="$PARAMETER_ROOT/$label.par"
readonly PROGRAM_LOG="$CAMPAIGN_ROOT/logs/$label.log"
readonly SCHEME="${reconstruction}_hlle"
mkdir -p "$PARAMETER_ROOT" "$CASE_ROOT" "$STATE_ROOT" "$CAMPAIGN_ROOT/logs"

if [[ -s "$STATE_ROOT/$label.complete" ]]; then
  printf 'SKIP %s: ya está completo.\n' "$label"
  exit 0
fi

cat >"$PARAMETER_FILE" <<EOF
# Campaña urgente para presentación: sensor activo.
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
nx = 300
ny = 1
nz = 150
r_min = 1.2
r_max = 40.0
theta_min = 0.0
theta_max = 1.0
phi_min = 0.0
phi_max = 2.0
logarithmic_r = true
reconstruction = $reconstruction
riemann_solver = hlle
tvd_limiter = mc
shock_sensor = true
final_time = 2000.0
cfl = 0.4
save_interval = 25.0
checkpoint_interval = 250.0
output_prefix = $label
output_folder = results/$label/session_000
vtk_mapping = physical
extract_gw = true
extract_mdot = true
ppi_diagnostics = true
diagnostic_stride = 500
extraction_stride = 500
mass_monitor_stride = 0
perturbation = pressure_noise
perturbation_seed = 3435
perturbation_time = 500.0
perturbation_amplitude = 0.01
perturbation_mode = 1.0
EOF

cd "$CAMPAIGN_ROOT"
"$EXECUTABLE" --check "$PARAMETER_FILE"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OMP_DYNAMIC=FALSE
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export GFORTRAN_UNBUFFERED_ALL=y

printf 'START %s case=%s job=%s task=%s host=%s threads=%s CPUs=%s\n' \
  "$(date --iso-8601=seconds)" "$label" "$SLURM_ARRAY_JOB_ID" "$TASK_ID" \
  "$(hostname)" "$OMP_NUM_THREADS" \
  "$(awk '/^Cpus_allowed_list:/ {print $2}' /proc/self/status)"

srun --ntasks=1 --cpus-per-task="$OMP_NUM_THREADS" --cpu-bind=cores \
  "$EXECUTABLE" "$PARAMETER_FILE" >"$PROGRAM_LOG" 2>&1

grep -q 'SIMULACIÓN FINALIZADA EXITOSAMENTE EN T' "$PROGRAM_LOG"
readonly DIAGNOSTICS="$CASE_ROOT/$SCHEME/global_diagnostics.dat"
readonly EVENT_FILE="$CASE_ROOT/$SCHEME/perturbation_events.dat"
test -s "$DIAGNOSTICS"
test -s "$EVENT_FILE"
last_time=$(awk '!/^#/ && NF {value=$1} END {print value}' "$DIAGNOSTICS")
event_count=$(awk '!/^#/ && NF {count++} END {print count+0}' "$EVENT_FILE")
awk -v value="$last_time" 'BEGIN {exit(value >= 1999.999999 && value <= 2000.000001 ? 0 : 1)}'
(( event_count == 1 ))

printf 'case=%s\njob=%s\ntask=%s\nthreads=%s\ncompleted=%s\n' \
  "$label" "$SLURM_ARRAY_JOB_ID" "$TASK_ID" "$OMP_NUM_THREADS" \
  "$(date --iso-8601=seconds)" >"$STATE_ROOT/$label.complete"
printf 'COMPLETE %s case=%s t=%s\n' "$(date --iso-8601=seconds)" "$label" "$last_time"
