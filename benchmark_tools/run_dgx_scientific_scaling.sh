#!/usr/bin/env bash
# The unchanged27-run scientific panel; one exclusive task at a time.
set -euo pipefail
root=$(realpath "${1:?Require project root}")
recipe_sha=${2:?Require frozen recipe manifest SHA256}
spec_sha=${3:?Require scientific specification SHA256}
index=${SLURM_ARRAY_TASK_ID:?Require array index}
[[ "$recipe_sha" =~ ^[a-f0-9]{64}$ && "$spec_sha" =~ ^[a-f0-9]{64}$ ]]
[[ "$index" =~ ^([0-9]|1[0-9]|2[0-6])$ ]]
recipe="$root/native_launcher_recipe_v2/benchmark_tools"
printf -v run 'run_%02d' "$index"
test ! -e "$root/scaling_native_v1/$run"
test ! -e "$root/scaling_native_v1/${run}_cache"
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
export PYTHONPYCACHEPREFIX="$root/scaling_native_v1/${run}_cache"
unset PYTHONPATH LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
exec "$root/envs/orthohmm/bin/python" -B "$recipe/launch_dgx_native_run.py" \
    --spec "$recipe/dgx_scientific_execution_20260917.json" --spec-sha256 "$spec_sha" \
    --recipe-manifest "$root/runtime_inventory_v1/native_launcher_recipe_v2.json" \
    --recipe-sha256 "$recipe_sha" --index "$index"
