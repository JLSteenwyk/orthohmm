#!/usr/bin/env bash
# Each invocation needs its own exclusive20CPU96GiB Slurm task allocation.
set -euo pipefail
root=$(realpath "${1:?Require project root}")
label=${2:?Require unique output label}
iterations=${3:?Require fixed iteration count}
interval=${4:?Require resource sampling interval}
[[ "$label" =~ ^collector_load_[a-zA-Z0-9_]+$ ]]
[[ "$iterations" =~ ^[1-9][0-9]*$ ]]
[[ "$interval" = 1 || "$interval" = 86400 ]]
test "$(hostname)" = spark-7ff0
test "${SLURM_CPUS_PER_TASK}" = 20
test ! -e "$root/$label"
cd "$root/collector_load_recipe_v2"
sha256sum --check <<'CHECKSUMS'
6da28892666f47053ae80ff66520897648f12a62413d739d918274c0f4c05cac  collector_load_fixture.c
80f97405b9361fa021b3857298ee7c4a3ee69e2f99d774cfba55a37afb0851ec  load
aaa0bda5d4a1966ffd712c92c6517555b57da3e710efe99ab84b18baa8b9ad14  measure_slurm_command.py
9e994aab177eff74c259416efdd509945708d1363271f63c996523447e346924  monitor_slurm_resources.py
3ef0f82132bad825509f8fade0e11c1b17d12be11d578327a7c593cbe3db30c9  command_host_monitor.py
d8b66001bdabf4c436aab8185f218dac019dd45b3075c71521513235a1c07e0a  observe_host_competition.py
4f78a03aab56c9e0771eba90ad15c4f91b3d5cb6289a9dd4b6407d444cec45d3  slurm_resource_snapshot.py
CHECKSUMS
export OMP_NUM_THREADS=20 OMP_DYNAMIC=FALSE OMP_PROC_BIND=TRUE OMP_PLACES=cores
unset PYTHONPATH GOMP_CPU_AFFINITY
exec "$root/envs/orthohmm/bin/python" measure_slurm_command.py \
    --output "$root/$label" --job-id "$SLURM_JOB_ID" --cpus 20 --memory-gib 96 \
    --timeout 600 --interval "$interval" --monitor-host --host-interval 30 \
    -- "$root/collector_load_recipe_v2/load" "$iterations" 20
