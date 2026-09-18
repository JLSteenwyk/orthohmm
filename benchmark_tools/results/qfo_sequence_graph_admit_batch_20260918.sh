#!/bin/bash
#SBATCH --job-name=qfo_sequence_graph_admit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=192G
#SBATCH --time=1-00:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_sequence_graph_admit_%j.log
set -euo pipefail
if [[ $# != 5 ]]; then
    printf '%s\n' 'Require PLAN PLAN_SHA VARIANT GRAPH_JOB REPORT_SHA' >&2
    exit 2
fi
PLAN=$1
PLAN_SHA=$2
VARIANT=$3
GRAPH_JOB=$4
REPORT_SHA=$5
if [[ "$PLAN" != /* || ! "$PLAN_SHA" =~ ^[a-f0-9]{64}$ || ! "$REPORT_SHA" =~ ^[a-f0-9]{64}$ || ! "$GRAPH_JOB" =~ ^[1-9][0-9]*$ ]]; then
    printf '%s\n' 'Invalid absolute plan path, SHA256 or graph job ID' >&2
    exit 2
fi
if [[ "$VARIANT" != all_hits && "$VARIANT" != top100 ]]; then
    printf '%s\n' 'Variant must be all_hits or top100' >&2
    exit 2
fi
if [[ -z "${SLURM_JOB_ID:-}" || "${SLURM_CPUS_PER_TASK:-}" != 2 || "${SLURM_MEM_PER_NODE:-}" != 196608 ]]; then
    printf '%s\n' 'Require scheduled 2CPU/192GiB admission task' >&2
    exit 2
fi
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR="$ROOT/benchmarks/work/publication_qfo_sequence_graph_admission_v1"
COMMIT=f1e21b09c28f270dc3ef2243bdcad86f212b58a0
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools orthohmm
test -f "$PLAN"
OUTPUT="$ROOT/benchmarks/work/qfo_sequence_graph_admission_${VARIANT}_20260918.json"
test ! -e "$OUTPUT"
test ! -L "$OUTPUT"
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/admit_qfo_sequence_graph.py" \
    --root "$ROOT" --plan "$PLAN" --plan-sha256 "$PLAN_SHA" --variant "$VARIANT" \
    --job "$GRAPH_JOB" --report-sha256 "$REPORT_SHA" --output "$OUTPUT"
