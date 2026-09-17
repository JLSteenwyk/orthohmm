#!/bin/bash
#SBATCH --job-name=resource_host_smoke
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/resource_host_smoke_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen executor required}
COMMIT=${2:?Frozen executor revision required}
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/measure_slurm_command.py" \
    --output "$ROOT/benchmarks/results/resource_host_command_smoke_v1" --job-id "$SLURM_JOB_ID" \
    --cpus 2 --memory-gib 1 --timeout 30 --interval 0.5 --monitor-host -- \
    /home/bizon/anaconda3/bin/python -c 'import subprocess,sys; a=bytearray(16*1024*1024); subprocess.run([sys.executable,"-c","import time; a=bytearray(32*1024*1024); sum(i*i for i in range(2000000)); time.sleep(3)"],check=True)'
