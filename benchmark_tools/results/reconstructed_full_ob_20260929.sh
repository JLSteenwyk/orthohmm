#!/bin/bash
#SBATCH --job-name=orthohmm_reconstructed_full_ob
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
  --plan "$ROOT/benchmarks/work/publication_reconstructed_full_ob_20260929/plan.json" \
  --sha256 cea0dd8efa459a01c57ca0450798005c5de3aeef15a6fda19c181acd35918bee
