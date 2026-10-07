#!/bin/bash
#SBATCH --job-name=ohmm_allocated_review10
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --no-requeue
set -euo pipefail
if [[ "${SLURM_CPUS_PER_TASK:-}" != 2 || ! "${SLURM_JOB_ID:-}" =~ ^[1-9][0-9]*$ || "${SLURM_JOB_ID}" == 23902 ]]; then
    printf '%s\n' 'Require a distinct scheduled two-CPU review job.' >&2
    exit 2
fi
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
cd "$ROOT"
PYTHON="$ROOT/benchmarks/work/native_factorial_review_py310_20261004/bin/python"
ENV=(env -u PYTHONHOME -u PYTHONPATH -u PYTHONUSERBASE
    -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT
    PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0
    PYTHONFAULTHANDLER=1 PYTHONUNBUFFERED=1
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1)
"${ENV[@]}" "$PYTHON" -B -c 'import json,time; from pathlib import Path; from benchmark_tools.run_native_factorial_cost import available_memory; from benchmark_tools.validate_native_factorial_outputs import require; raw=Path("/proc/meminfo").read_text(); n=available_memory(raw); print(json.dumps(dict(stage="allocated_native_terminal_review_capacity",observed_unix_ns=time.time_ns(),available_memory_bytes=n,minimum_available_memory_bytes=128*2**30,raw_meminfo=raw))); require(n >= 128*2**30, "Unsafe review capacity; retain this attempt without automatic retry")'
exec "${ENV[@]}" "$PYTHON" -B -X faulthandler \
    -m benchmark_tools.review_allocated_native_factorial_attempt \
    --request "$ROOT/benchmarks/work/native_factorial_launch_20261004/request_10_allocated_v1.json" \
    --request-sha256 1355ae3b73a9c8496133bd25a6cbf136d149598ccab680672e81f86482338910 \
    --output-directory "$ROOT/benchmarks/work/allocated_native_factorial_terminal_review_23902_v1"
