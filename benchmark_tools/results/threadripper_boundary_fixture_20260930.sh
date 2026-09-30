#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --exclusive
#SBATCH --time=00:05:00
#SBATCH --no-requeue
#SBATCH --job-name=oh-boundary-fixture
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/threadripper_boundary_fixture_20260930_%j.log
set -euo pipefail
cd /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
exec env -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/threadripper_private_controller_20260928/venv/bin/python \
  -B -m benchmark_tools.probe_threadripper_boundary \
  --output /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/threadripper_boundary_fixture_20260930 "$@"
