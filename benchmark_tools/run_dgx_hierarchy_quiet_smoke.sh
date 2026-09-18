#!/usr/bin/env bash
set -euo pipefail
root=/home/jlsteenwyk/projects/orthohmm-publication
recipe="$root/hierarchy_quiet_recipe_v1"
sha=${1:?Require recipe manifest SHA256}
index=${SLURM_ARRAY_TASK_ID:?Require array index}
[[ "$sha" =~ ^[a-f0-9]{64}$ && "$index" =~ ^[0-2]$ ]]
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
export PYTHONPYCACHEPREFIX="$root/hierarchy_quiet_smoke_v1/launcher_cache_$index"
test ! -e "$PYTHONPYCACHEPREFIX"
unset PYTHONPATH LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
exec "$root/envs/orthohmm/bin/python" -B "$recipe/benchmark_tools/run_dgx_hierarchy_quiet_smoke.py" \
  --spec "$recipe/benchmark_tools/dgx_launcher_smoke_spec_20260917.json" \
  --recipe "$root/hierarchy_quiet_recipe_v1.json" --recipe-sha "$sha" --index "$index"
