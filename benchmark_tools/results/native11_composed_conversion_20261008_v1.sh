#!/bin/bash
#SBATCH --job-name=ohmm_native11_convert
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --no-requeue
set -euo pipefail
if [[ "$#" != 5 || "${SLURM_CPUS_PER_TASK:-}" != 2 || ! "${SLURM_JOB_ID:-}" =~ ^[1-9][0-9]*$ ]]; then
    printf '%s\n' 'Require scheduled two-CPU conversion and request/review/source bindings.' >&2
    exit 2
fi
for digest in "$2" "$4" "$5"; do
    [[ "$digest" =~ ^[0-9a-f]{64}$ ]] || exit 2
done
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
cd "$ROOT"
exec env -u PYTHONHOME -u PYTHONPATH -u PYTHONUSERBASE \
    -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
    PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
    PYTHONFAULTHANDLER=1 PYTHONUNBUFFERED=1 \
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    "$ROOT/benchmarks/work/native_factorial_review_py310_20261004/bin/python" \
    -B -X faulthandler -m benchmark_tools.prepare_native11_composed_qfo_pairs \
    --request "$1" --request-sha256 "$2" --terminal-review "$3" \
    --terminal-review-sha256 "$4" --source-sha256 "$5"
