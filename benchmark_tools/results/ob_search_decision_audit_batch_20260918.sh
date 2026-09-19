#!/bin/bash
#SBATCH --job-name=ob_search_decision_audit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/ob_search_decision_audit_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen auditor required}
COMMIT=${2:?Exact auditor revision required}
JOB=${3:?Diagnostic job required}
[[ $# == 3 && "$COMMIT" =~ ^[0-9a-f]{40}$ && "$JOB" =~ ^[0-9]+$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --exit-code HEAD -- benchmark_tools
REPORT="$ROOT/benchmarks/results/ob_search_decisions_v1/report.json"
read -r SHA _ < <(sha256sum "$REPORT")
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/audit_ob_search_decisions.py" \
    --root "$ROOT" --report "$REPORT" --report-sha256 "$SHA" --job "$JOB" \
    --output "$ROOT/benchmarks/work/ob_search_decisions_audit_20260918.json"
