#!/usr/bin/env bash
# Sequential native-launcher smoke array, not the scientific scaling panel.
set -euo pipefail
root=$(realpath "${1:?Require project root}")
recipe_sha=${2:?Require frozen recipe manifest SHA256}
spec_sha=${3:?Require authorized smoke specification SHA256}
index=${SLURM_ARRAY_TASK_ID:?Require array index}
[[ "$recipe_sha" =~ ^[a-f0-9]{64}$ && "$spec_sha" =~ ^[a-f0-9]{64}$ ]]
[[ "$index" =~ ^[0-2]$ ]]
recipe="$root/native_launcher_recipe_v1/benchmark_tools"
printf -v run 'run_%02d' "$index"
test ! -e "$root/launcher_smoke_v1/$run"
test ! -e "$root/launcher_smoke_v1/${run}_cache"
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
export PYTHONPYCACHEPREFIX="$root/launcher_smoke_v1/${run}_cache"
unset PYTHONPATH LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
exec "$root/envs/orthohmm/bin/python" -B "$recipe/launch_dgx_native_run.py" \
    --spec "$recipe/dgx_launcher_smoke_spec_20260917.json" --spec-sha256 "$spec_sha" \
    --recipe-manifest "$root/runtime_inventory_v1/native_launcher_recipe_v1.json" \
    --recipe-sha256 "$recipe_sha" --index "$index"
