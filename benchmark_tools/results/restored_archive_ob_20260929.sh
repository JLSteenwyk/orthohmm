#!/bin/bash
#SBATCH --job-name=orthohmm_restored_archive_ob
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --nodelist=bizon
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
exec /tmp/orthohmm-base-reconstruction-20260927/python-runtime/bin/python -I -S -B \
  "$ROOT/benchmark_tools/run_integrated_full_job.py" \
  --plan "$ROOT/benchmarks/work/publication_restored_archive_ob_20260929/plan.json" \
  --sha256 e80067f7a5947bc357831e0d65c7d766d09899f004cfc6a31d1274a8dfe15448
