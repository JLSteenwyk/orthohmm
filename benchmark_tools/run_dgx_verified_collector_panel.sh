#!/usr/bin/env bash
# Frozen counterbalanced engineering panel; never scientific inference.
set -euo pipefail
root=$(realpath "${1:?Require project root}")
recipe_sha=${2:?Require prospective recipe-tree manifest SHA256}
[[ "$recipe_sha" =~ ^[a-f0-9]{64}$ ]]
index=${SLURM_ARRAY_TASK_ID:?Require array index}
case "$index" in
    0|3|4) mode=sparse; interval=86400 ;;
    1|2|5) mode=sampled; interval=1 ;;
    *) exit 2 ;;
esac
recipe="$root/verified_collector_recipe_v2"
output="$root/verified_collector_panel_v2_${index}_${mode}"
test "$(hostname)" = spark-7ff0
test "${SLURM_CPUS_PER_TASK}" = 20
test ! -e "$output"
test ! -e "${output}_cache"
cd "$recipe"
sha256sum --check <<'CHECKSUMS'
36243fcdfe6f49e73dbe1dd1b627bf4a330f19f3cd52dcf94a9618ddf0415cf8  run_verified_slurm_measurement.py
ab0008651e75ebadc1be3c87503f09bffd9e5be87f2e6b0c3eec49e552d9bd4c  snapshot_runtime_trees.py
CHECKSUMS
export PYTHONPYCACHEPREFIX="${output}_cache" PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=20 OMP_DYNAMIC=FALSE OMP_PROC_BIND=TRUE OMP_PLACES=cores
unset PYTHONPATH GOMP_CPU_AFFINITY LD_PRELOAD LD_LIBRARY_PATH LD_AUDIT
exec "$root/envs/orthohmm/bin/python" -B "$recipe/run_verified_slurm_measurement.py" \
    --runtime "$root/runtime_inventory_v1/trees.json" 2f38fc57683e51a7b6768b4709db16293590983640c41cfe12a56ed6036750ae \
    --runtime "$root/runtime_inventory_v1/system_trees.json" 4083d15c0ffc40013756c4588763aa8af24ca39165622665ee5074662e4b46ce \
    --runtime "$root/runtime_inventory_v1/recipe_trees_v2.json" "$recipe_sha" \
    --collector "$root/collector_load_recipe_v2" --cwd "$root/collector_load_recipe_v2" \
    --output "$output" --job-id "$SLURM_JOB_ID" --cpus 20 --memory-gib 96 --timeout 600 --interval "$interval" \
    -- /usr/bin/time -q -f $'elapsed_seconds\t%e\nuser_seconds\t%U\nsystem_seconds\t%S\nmax_process_rss_kib\t%M\nexit_status\t%x' \
    -o "$output/native.time.tsv" -- "$root/collector_load_recipe_v2/load" 8000000000 20
