#!/usr/bin/env bash
# Submit as a six-task array with concurrency one; no scientific inference.
set -euo pipefail
root=/home/jlsteenwyk/projects/orthohmm-publication
cd "$root/collector_load_recipe_v2"
sha256sum --check <<'CHECKSUMS'
a323dacc0b0df4443b3dbdce3b1a786f6e4293ecac153de291e868d3a5aa8fca  run_dgx_collector_load.sh
CHECKSUMS
index=${SLURM_ARRAY_TASK_ID:?Require array identity}
case "$index" in
    0|3|4) mode=sparse; interval=86400 ;;
    1|2|5) mode=sampled; interval=1 ;;
    *) exit 2 ;;
esac
exec bash "$root/collector_load_recipe_v2/run_dgx_collector_load.sh" "$root" \
    "collector_load_panel_v1_${index}_${mode}" 8000000000 "$interval"
