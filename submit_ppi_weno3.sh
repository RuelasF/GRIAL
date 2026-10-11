#!/usr/bin/env bash
set -Eeuo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
walltime=${WALLTIME:-3-00:00:00}
cpus_per_case=${CPUS_PER_CASE:-8}
mail_user=${MAIL_USER:-francisco.ruelas9565@alumnos.udg.mx}
ppi_nx=${PPI_NX:-400}
ppi_nz=${PPI_NZ:-200}

if ! command -v sbatch >/dev/null 2>&1; then
  printf 'ERROR: sbatch no esta disponible. Ejecute este script en el cluster.\n' >&2
  exit 1
fi
if [[ ! "$cpus_per_case" =~ ^[1-9][0-9]*$ ]]; then
  printf 'ERROR: CPUS_PER_CASE debe ser un entero positivo.\n' >&2
  exit 1
fi
if [[ ! "$ppi_nx" =~ ^[1-9][0-9]*$ || ! "$ppi_nz" =~ ^[1-9][0-9]*$ ]]; then
  printf 'ERROR: PPI_NX y PPI_NZ deben ser enteros positivos.\n' >&2
  exit 1
fi

mkdir -p "$script_dir/slurm_logs"

if [[ ! -x "$script_dir/grial" ]]; then
  printf 'ERROR: no existe el ejecutable %s/grial. Compilelo antes de enviar.\n' \
    "$script_dir" >&2
  exit 1
fi

declare -a site_options=()
[[ -z "${ACCOUNT:-}" ]] || site_options+=(--account="$ACCOUNT")
[[ -z "${PARTITION:-}" ]] || site_options+=(--partition="$PARTITION")
[[ -z "${RESERVATION:-}" ]] || site_options+=(--reservation="$RESERVATION")

batch_one_submission=$(sbatch --parsable \
  --job-name=ppi-w3-b1 \
  --array=0-2 \
  --nodes=1 --ntasks=1 --cpus-per-task="$cpus_per_case" \
  --time="$walltime" \
  --chdir="$script_dir" \
  --export="ALL,GRHD_REPO_ROOT=$script_dir,PPI_NX=$ppi_nx,PPI_NZ=$ppi_nz" \
  --output="$script_dir/slurm_logs/weno3_%A_%a.log" \
  --mail-user="$mail_user" \
  --mail-type=BEGIN,END,FAIL,ARRAY_TASKS \
  "${site_options[@]}" \
  "$script_dir/ppi_weno3_case.sbatch")
batch_one=${batch_one_submission%%;*}

batch_two_submission=$(sbatch --parsable \
  --dependency="afterok:$batch_one" \
  --job-name=ppi-w3-b2 \
  --array=3-5 \
  --nodes=1 --ntasks=1 --cpus-per-task="$cpus_per_case" \
  --time="$walltime" \
  --chdir="$script_dir" \
  --export="ALL,GRHD_REPO_ROOT=$script_dir,PPI_NX=$ppi_nx,PPI_NZ=$ppi_nz" \
  --output="$script_dir/slurm_logs/weno3_%A_%a.log" \
  --mail-user="$mail_user" \
  --mail-type=BEGIN,END,FAIL,ARRAY_TASKS \
  "${site_options[@]}" \
  "$script_dir/ppi_weno3_case.sbatch")
batch_two=${batch_two_submission%%;*}

printf 'Ejecutable:         %s/grial (precompilado)\n' "$script_dir"
printf 'Malla:              %sx1x%s\n' "$ppi_nx" "$ppi_nz"
printf 'Primera tanda 0-2: %s\n' "$batch_one"
printf 'Segunda tanda 3-5: %s (depende de %s)\n' "$batch_two" "$batch_one"
printf 'Seguimiento: squeue -j %s,%s\n' "$batch_one" "$batch_two"
printf 'Logs Slurm:  %s/slurm_logs/\n' "$script_dir"
