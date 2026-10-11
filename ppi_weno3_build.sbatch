#!/usr/bin/env bash
#SBATCH --job-name=grhd-build
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --output=build_%j.log

# La asignacion de CPU, tiempo, cuenta y particion se hace desde submit_ppi_weno3.sh.

set -Eeuo pipefail

repo_root=${GRHD_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}
cd "$repo_root"

[[ -f Makefile ]] || {
  printf 'ERROR: no se encontro Makefile en GRHD_REPO_ROOT=%s.\n' "$repo_root" >&2
  exit 1
}

command -v gfortran >/dev/null 2>&1 || {
  printf 'ERROR: gfortran no esta en PATH. Cargue el modulo del compilador antes de enviar.\n' >&2
  exit 1
}

printf 'Nodo de compilacion: %s\n' "$(hostname)"
if command -v git >/dev/null 2>&1 && [[ -d .git ]]; then
  printf 'Commit: %s\n' "$(git rev-parse HEAD)"
else
  printf 'Commit: no disponible (git no esta instalado en el nodo)\n'
fi
gfortran --version | head -n 1

make clean
make -j "${SLURM_CPUS_PER_TASK:-1}" grial
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_DYNAMIC=FALSE
export OMP_PLACES=cores
export OMP_PROC_BIND=close
make test

printf 'Compilacion y pruebas OpenMP terminadas: %s/grial\n' "$repo_root"
