#!/bin/bash
#SBATCH --job-name=resource_smoke
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/resource_command_smoke_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$ROOT/benchmarks/work/publication_resource_command_v1/benchmark_tools/measure_slurm_command.py" --output "$ROOT/benchmarks/results/resource_command_smoke_v1" --job-id "$SLURM_JOB_ID" --cpus 2 --memory-gib 1 --timeout 30 --interval 0.2 -- /home/bizon/anaconda3/bin/python -c 'import subprocess,sys; a=bytearray(16*1024*1024); subprocess.run([sys.executable,"-c","import time; a=bytearray(32*1024*1024); sum(i*i for i in range(2000000)); time.sleep(2)"],check=True)'
