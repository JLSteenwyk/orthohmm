#!/usr/bin/env bash
# Engineering smoke only: runtime checks bracket the unchanged native collector.
set -euo pipefail
root=$(realpath "${1:?Require project root}")
recipe=$(realpath "${2:?Require recipe directory}")
test "$(hostname)" = spark-7ff0
test "${SLURM_CPUS_PER_TASK}" = 20
test ! -e "$root/verified_collector_smoke_v1"
test ! -e "$root/verified_collector_smoke_v1_cache"
cd "$recipe"
sha256sum --check <<'CHECKSUMS'
36243fcdfe6f49e73dbe1dd1b627bf4a330f19f3cd52dcf94a9618ddf0415cf8  run_verified_slurm_measurement.py
ab0008651e75ebadc1be3c87503f09bffd9e5be87f2e6b0c3eec49e552d9bd4c  snapshot_runtime_trees.py
CHECKSUMS
export PYTHONPYCACHEPREFIX="$root/verified_collector_smoke_v1_cache" PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=20 OMP_DYNAMIC=FALSE OMP_PROC_BIND=TRUE OMP_PLACES=cores
unset PYTHONPATH GOMP_CPU_AFFINITY LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
exec "$root/envs/orthohmm/bin/python" -B "$recipe/run_verified_slurm_measurement.py" \
    --runtime "$root/runtime_inventory_v1/trees.json" 2f38fc57683e51a7b6768b4709db16293590983640c41cfe12a56ed6036750ae \
    --runtime "$root/runtime_inventory_v1/system_trees.json" 4083d15c0ffc40013756c4588763aa8af24ca39165622665ee5074662e4b46ce \
    --collector "$root/collector_load_recipe_v2" --cwd "$root/collector_load_recipe_v2" \
    --output "$root/verified_collector_smoke_v1" --job-id "$SLURM_JOB_ID" \
    --cpus 20 --memory-gib 96 --timeout 300 --interval 1 \
    -- /usr/bin/time -q -f $'elapsed_seconds\t%e\nuser_seconds\t%U\nsystem_seconds\t%S\nmax_process_rss_kib\t%M\nexit_status\t%x' \
    -o "$root/verified_collector_smoke_v1/native.time.tsv" -- \
    "$root/collector_load_recipe_v2/load" 1000000000 20
