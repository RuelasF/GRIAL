#!/usr/bin/env bash
#SBATCH --job-name=jet-omp-scale
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:30:00
#SBATCH --output=jet_scaling_%j.log

# Una sola asignacion de 16 CPU ejecuta 1, 2, 4, 8 y 16 hilos secuencialmente.

set -Eeuo pipefail
umask 027

repo_root=${GRHD_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}
parameter_file="$repo_root/benchmarks/jet/jet_openmp_scaling.par"
executable="$repo_root/grial"
result_root="$repo_root/jet_openmp_scaling/job_${SLURM_JOB_ID:-manual}"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OMP_DYNAMIC=FALSE
export OMP_PLACES=cores
export OMP_PROC_BIND=close

cd "$repo_root"
[[ -f Makefile ]] || {
  printf 'ERROR: no se encontro Makefile en %s.\n' "$repo_root" >&2
  exit 1
}
[[ -f "$parameter_file" ]] || {
  printf 'ERROR: no se encontro %s.\n' "$parameter_file" >&2
  exit 1
}
command -v gfortran >/dev/null 2>&1 || {
  printf 'ERROR: gfortran no esta disponible en el nodo.\n' >&2
  exit 1
}
command -v python3 >/dev/null 2>&1 || {
  printf 'ERROR: python3 no esta disponible para comparar los VTK.\n' >&2
  exit 1
}

mkdir -p "$result_root"
printf 'Nodo: %s\n' "$(hostname)"
printf 'Directorio del repositorio: %s\n' "$repo_root"
printf 'CPUs de la asignacion: %s\n' "${SLURM_CPUS_PER_TASK:-16}"
printf 'OpenMP recibe inicialmente OMP_NUM_THREADS=%s\n' "$OMP_NUM_THREADS"
printf 'OMP_PLACES=%s; OMP_PROC_BIND=%s\n' "$OMP_PLACES" "$OMP_PROC_BIND"
printf 'CPUs permitidas al batch: %s\n' \
  "$(awk '/^Cpus_allowed_list:/ {print $2}' /proc/self/status)"
gfortran --version | head -n 1

make -j "${SLURM_CPUS_PER_TASK:-16}" grial
"$executable" --check "$parameter_file"

thread_counts=(1 2 4 8 16)
for threads in "${thread_counts[@]}"; do
  run_root="$result_root/threads_$threads"
  mkdir -p "$run_root"
  printf '\n===== JET CON %s HILOS OPENMP =====\n' "$threads"
  start_ns=$(date +%s%N)
  (
    cd "$run_root"
    srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task="$threads" --cpu-bind=cores \
      env OMP_NUM_THREADS="$threads" OMP_DYNAMIC=FALSE OMP_PLACES=cores \
      OMP_PROC_BIND=close OMP_DISPLAY_ENV=TRUE OMP_DISPLAY_AFFINITY=TRUE \
      GFORTRAN_UNBUFFERED_ALL=y \
      "$executable" "$parameter_file" >run.log 2>&1
  )
  finish_ns=$(date +%s%N)
  awk -v start="$start_ns" -v finish="$finish_ns" \
    'BEGIN {printf "%.9f\n", (finish-start)/1000000000}' \
    >"$run_root/wall_seconds.txt"

  final_vtk="$run_root/jet_scaling_output/weno5_hlle/JET_OMP_weno5_hlle_0001.vtk"
  [[ -s "$final_vtk" ]] || {
    printf 'ERROR: la corrida con %s hilos no genero el VTK final.\n' "$threads" >&2
    exit 1
  }
  grep -q 'FINALIZADA EXITOSAMENTE' "$run_root/run.log" || {
    printf 'ERROR: la corrida con %s hilos no informo terminacion correcta.\n' "$threads" >&2
    exit 1
  }
  printf 'Tiempo con %s hilos: %s s\n' "$threads" "$(<"$run_root/wall_seconds.txt")"
done

reference="$result_root/threads_1/jet_scaling_output/weno5_hlle/JET_OMP_weno5_hlle_0001.vtk"
comparison_status=0
for threads in "${thread_counts[@]:1}"; do
  candidate="$result_root/threads_$threads/jet_scaling_output/weno5_hlle/JET_OMP_weno5_hlle_0001.vtk"
  python3 "$repo_root/benchmarks/jet/compare_vtk.py" "$reference" "$candidate" \
    --relative-tolerance 1e-10 --absolute-tolerance 1e-12 \
    >"$result_root/compare_1_vs_${threads}.txt" || comparison_status=$?
done

baseline=$(<"$result_root/threads_1/wall_seconds.txt")
{
  printf 'threads\twall_seconds\tspeedup\tefficiency\n'
  for threads in "${thread_counts[@]}"; do
    value=$(<"$result_root/threads_$threads/wall_seconds.txt")
    awk -v base="$baseline" -v value="$value" -v count="$threads" \
      'BEGIN {speed=base/value; printf "%d\t%.9f\t%.6f\t%.6f\n", count, value, speed, speed/count}'
  done
} | tee "$result_root/timing.tsv"

printf 'Resultados: %s\n' "$result_root"
if (( comparison_status != 0 )); then
  printf 'ERROR: al menos una comparacion numerica excedio la tolerancia.\n' >&2
  exit "$comparison_status"
fi
