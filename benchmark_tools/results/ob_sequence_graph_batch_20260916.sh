#!/bin/bash
#SBATCH --job-name=ob_seq_graph
#SBATCH --array=0-1%1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/ob_seq_graph_%A_%a.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
case "$SLURM_ARRAY_TASK_ID" in
  0) VARIANT=all_hits ;;
  1) VARIANT=top100 ;;
  *) exit 2 ;;
esac
exec /home/bizon/anaconda3/bin/python "$ROOT/benchmarks/work/publication_ob_sequence_graph_v1/benchmark_tools/run_sequence_graph_control.py" --root "$ROOT" --variant "$VARIANT"
