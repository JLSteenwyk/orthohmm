#!/usr/bin/env bash
# Read-only engineering probe; no native inference or synthetic CPU load.
set -euo pipefail
recipe=$(realpath "${1:?Require probe directory}")
test "$(hostname)" = spark-7ff0
test "${SLURM_CPUS_PER_TASK}" = 1
test "${SLURM_MEM_PER_NODE}" = 1024
test -n "${SLURM_JOB_ID}"
cd "$recipe"
sha256sum --check <<'CHECKSUMS'
3d3c193305085e17aa10e05be842d814f66d020d7dde285c8d9288ead2ce3f3c  probe_host_counters.py
CHECKSUMS
test ! -e "slurm_${SLURM_JOB_ID}.json"
printf 'job_id=%s\nnode=%s\ncpus=%s\nmemory_mib=%s\n' \
    "$SLURM_JOB_ID" "$(hostname)" "$SLURM_CPUS_PER_TASK" "$SLURM_MEM_PER_NODE"
exec /usr/bin/python3 -I -B "$recipe/probe_host_counters.py" \
    --interval 3 --output "$recipe/slurm_${SLURM_JOB_ID}.json"
