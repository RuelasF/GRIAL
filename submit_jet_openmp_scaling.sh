#!/usr/bin/env bash
set -Eeuo pipefail

repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
walltime=${WALLTIME:-00:30:00}

command -v sbatch >/dev/null 2>&1 || {
  printf 'ERROR: sbatch no esta disponible. Ejecute este script en el cluster.\n' >&2
  exit 1
}
mkdir -p "$repo_root/slurm_logs"

declare -a site_options=()
[[ -z "${ACCOUNT:-}" ]] || site_options+=(--account="$ACCOUNT")
[[ -z "${PARTITION:-}" ]] || site_options+=(--partition="$PARTITION")
[[ -z "${RESERVATION:-}" ]] || site_options+=(--reservation="$RESERVATION")

submission=$(sbatch --parsable \
  --job-name=jet-omp-scale \
  --nodes=1 --ntasks=1 --cpus-per-task=16 \
  --time="$walltime" \
  --chdir="$repo_root" \
  --export="ALL,GRHD_REPO_ROOT=$repo_root" \
  --output="$repo_root/slurm_logs/jet_scaling_%j.log" \
  "${site_options[@]}" \
  "$repo_root/jet_openmp_scaling.sbatch")
job_id=${submission%%;*}

printf 'Trabajo enviado: %s\n' "$job_id"
printf 'Estado: squeue -j %s\n' "$job_id"
printf 'Log: tail -f %s/slurm_logs/jet_scaling_%s.log\n' "$repo_root" "$job_id"
