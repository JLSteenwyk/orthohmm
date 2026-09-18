#!/usr/bin/env bash
set -euo pipefail
recipe=$(realpath "${1:?Require recipe directory}")
test "$(hostname)" = spark-7ff0
test "${SLURM_CPUS_PER_TASK}" = 2
test "${SLURM_MEM_PER_NODE}" = 2048
cd "$recipe/benchmark_tools"
sha256sum --check <<'CHECKSUMS'
e2ae95faad1630f0c055d0b5a152b87991e0c8a529ec58d7a43e072b25ee55e1  probe_dgx_cpu_hierarchy.py
3d3c193305085e17aa10e05be842d814f66d020d7dde285c8d9288ead2ce3f3c  probe_host_counters.py
c00eae1c44c86f0ee02c72e4634cf4e5cc70e2656a8ebfc0f9070769f6f7f97b  probe_dgx_step_separation.py
ef1040d6b588fc0ffaa2e554af9dcd9e8045b3c03874f9ad9d3dc1a3af627b70  DGX_CPU_HIERARCHY_PROTOCOL_20260918.md
CHECKSUMS
exec /usr/bin/python3 -I -B "$recipe/benchmark_tools/probe_dgx_cpu_hierarchy.py" --output "$recipe/job_${SLURM_JOB_ID}"
