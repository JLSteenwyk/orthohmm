#!/usr/bin/env bash
set -euo pipefail
recipe=$(realpath "${1:?Require recipe directory}")
test "$(hostname)" = spark-7ff0
test "${SLURM_CPUS_PER_TASK}" = 2
test "${SLURM_MEM_PER_NODE}" = 2048
cd "$recipe"
sha256sum --check <<'CHECKSUMS'
3d3c193305085e17aa10e05be842d814f66d020d7dde285c8d9288ead2ce3f3c  probe_host_counters.py
c00eae1c44c86f0ee02c72e4634cf4e5cc70e2656a8ebfc0f9070769f6f7f97b  probe_dgx_step_separation.py
b05e3750ba30692c728bbbeb166699e54bb5c949e6771786a61750942b6027eb  probe_dgx_compute_counters.py
56971dd3b0b9771ef2f036dc374267a8e9cbfbc7f920e0b178b4702b0588d6f4  DGX_COMPUTE_COUNTER_PROTOCOL_20260918.md
CHECKSUMS
exec /usr/bin/python3 -I -B "$recipe/probe_dgx_compute_counters.py" --output "$recipe/job_${SLURM_JOB_ID}"
