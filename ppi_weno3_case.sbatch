#!/usr/bin/env bash
#SBATCH --job-name=ppi-w3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00
#SBATCH --output=weno3_%A_%a.log

# Una tarea del arreglo Slurm ejecuta exactamente un caso WENO3.

set -Eeuo pipefail
umask 027

repo_root=${GRHD_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}
nx=${PPI_NX:-400}
nz=${PPI_NZ:-200}
campaign_root=${PPI_WENO3_ROOT:-"$repo_root/ppi_weno3_${nx}x1x${nz}_t5000"}
executable=${GRHD_EXECUTABLE:-"$repo_root/grial"}
task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID no esta definido}

if [[ ! "$nx" =~ ^[1-9][0-9]*$ || ! "$nz" =~ ^[1-9][0-9]*$ ]]; then
  printf 'ERROR: PPI_NX y PPI_NZ deben ser enteros positivos.\n' >&2
  exit 2
fi

cases=(
  ef_a00_weno3
  ks_a00_weno3
  ks_a02_weno3
  ks_a04_weno3
  ks_a06_weno3
  ks_a09_weno3
)

if (( task_id < 0 || task_id >= ${#cases[@]} )); then
  printf 'ERROR: indice de arreglo fuera de rango: %s\n' "$task_id" >&2
  exit 2
fi

label=${cases[$task_id]}
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
  ks_a02_weno3) metric=kerr_schild; geometry=spheroidal; spin=0.2 ;;
  ks_a04_weno3) metric=kerr_schild; geometry=spheroidal; spin=0.4 ;;
  ks_a06_weno3) metric=kerr_schild; geometry=spheroidal; spin=0.6 ;;
  ks_a09_weno3) metric=kerr_schild; geometry=spheroidal; spin=0.9 ;;
esac

parameter_root="$campaign_root/parameters"
result_root="$campaign_root/results/$label"
state_root="$campaign_root/state/$label"
mkdir -p "$parameter_root" "$result_root" "$state_root"

if [[ -s "$state_root/complete" ]]; then
  printf 'SKIP %s: el caso ya esta marcado como completo.\n' "$label"
  exit 0
fi

[[ -x "$executable" ]] || {
  printf 'ERROR: no existe el ejecutable %s. Revise el trabajo de compilacion.\n' "$executable" >&2
  exit 1
}

latest_checkpoint=$(find "$result_root" -type f -name 'checkpoint_*.rst' \
  -printf '%T@ %p\n' 2>/dev/null | sort -n | tail -n 1 | awk '{print $2}')

if [[ -n "$latest_checkpoint" ]]; then
  previous_session=$(find "$result_root" -mindepth 1 -maxdepth 1 -type d \
    -name 'session_[0-9][0-9][0-9]' -printf '%f\n' | sort | tail -n 1)
  previous_number=${previous_session#session_}
  session=$(printf 'session_%03d' "$((10#$previous_number + 1))")
  parameter_file="$parameter_root/${label}_${session}.par"
  cat >"$parameter_file" <<EOF
restart = true
restart_file = $latest_checkpoint
final_time = 5000.0
cfl = 0.4
save_interval = 50.0
checkpoint_interval = 250.0
output_prefix = $label
output_folder = results/$label/$session
vtk_mapping = physical
shock_sensor = false
extract_gw = true
extract_mdot = true
ppi_diagnostics = true
diagnostic_stride = 1000
extraction_stride = 1000
mass_monitor_stride = 0
perturbation = pressure_noise
perturbation_seed = 3435
perturbation_time = 1000.0
perturbation_amplitude = 0.01
perturbation_mode = 1.0
EOF
else
  session=session_000
  parameter_file="$parameter_root/${label}.par"
  cat >"$parameter_file" <<EOF
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
nx = $nx
ny = 1
nz = $nz
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
final_time = 5000.0
cfl = 0.4
save_interval = 50.0
checkpoint_interval = 250.0
output_prefix = $label
output_folder = results/$label/$session
vtk_mapping = physical
extract_gw = true
extract_mdot = true
ppi_diagnostics = true
diagnostic_stride = 1000
extraction_stride = 1000
mass_monitor_stride = 0
perturbation = pressure_noise
perturbation_seed = 3435
perturbation_time = 1000.0
perturbation_amplitude = 0.01
perturbation_mode = 1.0
EOF
fi

cd "$campaign_root"
"$executable" --check "$parameter_file"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OMP_DYNAMIC=FALSE
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export GFORTRAN_UNBUFFERED_ALL=y

printf 'Caso: %s\n' "$label"
printf 'Sesion: %s\n' "$session"
printf 'Nodo: %s\n' "$(hostname)"
printf 'Hilos OpenMP: %s\n' "$OMP_NUM_THREADS"
printf 'CPUs permitidas: %s\n' "$(awk '/^Cpus_allowed_list:/ {print $2}' /proc/self/status)"

srun --ntasks=1 --cpus-per-task="$OMP_NUM_THREADS" --cpu-bind=cores \
  "$executable" "$parameter_file"

diagnostics="$campaign_root/results/$label/$session/weno3_hlle/global_diagnostics.dat"
if [[ ! -s "$diagnostics" ]]; then
  printf 'ERROR: no se genero el diagnostico final de %s.\n' "$label" >&2
  exit 1
fi
last_time=$(awk '!/^#/ && NF {value=$1} END {print value}' "$diagnostics")
awk -v value="$last_time" 'BEGIN {exit(value >= 4999.999999 && value <= 5000.000001 ? 0 : 1)}' || {
  printf 'ERROR: %s termino con diagnostico t=%s, no t=5000.\n' "$label" "$last_time" >&2
  exit 1
}

printf 'case=%s\nsession=%s\njob_id=%s\nthreads=%s\ncompleted=%s\n' \
  "$label" "$session" "${SLURM_JOB_ID:-unknown}" "$OMP_NUM_THREADS" \
  "$(date --iso-8601=seconds)" >"$state_root/complete"
