#!/usr/bin/env bash
set -euo pipefail
root=/home/jlsteenwyk/projects/orthohmm-publication
recipe="$root/frontier_overhead_recipe_v1"
sha=${1:?Require recipe manifest SHA256}
auth=${2:?Require authorization SHA256}
index=${SLURM_ARRAY_TASK_ID:?Require array index}
[[ "$sha" =~ ^[a-f0-9]{64}$ && "$auth" =~ ^[a-f0-9]{64}$ && "$index" =~ ^([0-9]|1[0-7])$ ]]
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
export PYTHONPYCACHEPREFIX="$root/frontier_overhead_v1/launcher_cache_$index"
test ! -e "$PYTHONPYCACHEPREFIX"
unset PYTHONPATH LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
exec "$root/envs/orthohmm/bin/python" -B "$recipe/benchmark_tools/run_dgx_frontier_overhead.py" \
  --plan "$recipe/benchmark_tools/dgx_frontier_overhead_plan_20260918.json" \
  --authorization "$root/frontier_overhead_authorization_v1.json" --authorization-sha "$auth" \
  --recipe "$root/frontier_overhead_recipe_v1.json" --recipe-sha "$sha" --index "$index"
