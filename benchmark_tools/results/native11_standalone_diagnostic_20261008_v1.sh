#!/bin/bash
#SBATCH --job-name=ohmm_native11_diagnostic
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --no-requeue
set -euo pipefail
if [[ "$#" != 0 || "${SLURM_CPUS_PER_TASK:-}" != 2 || ! "${SLURM_JOB_ID:-}" =~ ^[1-9][0-9]*$ || "${SLURM_JOB_ID}" == 23985 || "${SLURM_JOB_ID}" == 23986 ]]; then
    printf '%s\n' 'Require a distinct scheduled two-CPU diagnostic job and no arguments.' >&2
    exit 2
fi
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
DEST="$ROOT/benchmarks/work/native11_standalone_diagnostic_20261008_v1"
cd "$ROOT"
test ! -e "$DEST"
test ! -L "$DEST"
PYTHON="$ROOT/benchmarks/work/native_factorial_review_py310_20261004/bin/python"
ENV=(env -u PYTHONHOME -u PYTHONPATH -u PYTHONUSERBASE
    -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT
    PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0
    PYTHONFAULTHANDLER=1 PYTHONUNBUFFERED=1
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1)
"${ENV[@]}" "$PYTHON" -B -c 'import json,time; from pathlib import Path; from benchmark_tools.run_native_factorial_cost import available_memory; from benchmark_tools.validate_native_factorial_outputs import require; raw=Path("/proc/meminfo").read_text(); n=available_memory(raw); print(json.dumps(dict(stage="native11_standalone_diagnostic_capacity",observed_unix_ns=time.time_ns(),available_memory_bytes=n,minimum_available_memory_bytes=32*2**30,raw_meminfo=raw))); require(n >= 32*2**30, "Unsafe diagnostic capacity; retain this attempt without automatic retry")'
mkdir "$DEST"
exec /usr/bin/time -v -o "$DEST/time.txt" "${ENV[@]}" "$PYTHON" -B -X faulthandler \
    -m benchmark_tools.validate_allocated_native_factorial_outputs \
    --request "$ROOT/benchmarks/work/native_factorial_launch_20261004/request_11_allocated_v1.json" \
    --request-sha256 7bf63b80bd5932b9edbd1b2c5ff3fb77f5557e4f6c64e077045d6a50c8d366a1 \
    --output "$DEST/outputs.json"
