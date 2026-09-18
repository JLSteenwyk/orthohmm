#!/usr/bin/env bash
set -euo pipefail
recipe=$(realpath "${1:?Require recipe directory}")
test "$(hostname)" = spark-7ff0
test "${SLURM_CPUS_PER_TASK}" = 2
test "${SLURM_MEM_PER_NODE}" = 2048
cd "$recipe/benchmark_tools"
sha256sum --check <<'CHECKSUMS'
3d3c193305085e17aa10e05be842d814f66d020d7dde285c8d9288ead2ce3f3c  probe_host_counters.py
c00eae1c44c86f0ee02c72e4634cf4e5cc70e2656a8ebfc0f9070769f6f7f97b  probe_dgx_step_separation.py
9175b0033b079d50cb086d64397a4996850d4f79dac317cd15ea83d96f4b8075  screen_bracketed_cpu.py
b535124efb386265555e1ea946e4edc16ad6d8ff1c3d467955576875d944e366  probe_interval_cpu.py
7798aab17925701cd4f315778e700c804dee9057b4a2dca0754c8011a9d2361f  DGX_INTERVAL_CONTROLS_PROTOCOL_20260918.md
CHECKSUMS
exec /usr/bin/python3 -I -B "$recipe/benchmark_tools/probe_interval_cpu.py" --output "$recipe/job_${SLURM_JOB_ID}"
