# Parameter Count-Audit Integration Check

Submitted job 21953 on 2026-09-19 at 09:00 EDT, with 2 CPUs, 32G RAM,
one-hour limit, node bizon and no requeue. Confirmed RUNNING at the initial
controller check. This is a control-only integration check, not a parameter
result or a complete robustness analysis.

Frozen executor: `benchmarks/work/publication_qfo_parameter_swiss_audit_v1`,
commit `99a09000014613c0906f8ed7331273b3c9759864` (pushed).
Input inventory SHA-256:
`25fbd435e7afb1dddd410433dc6df343c3c4586f3fa81a45980c39a19cdefa8c`.
The corrected full-pipeline control score admission and the historical
reference-universe audit are separately hard-pinned by the auditor. The
historical file supplies reference labels/members only, not prediction counts.

Slurm wrap command, run from the original repository root:

```bash
env -u PYTHONPATH -u PYTHONHOME -u LD_PRELOAD -u LD_LIBRARY_PATH \
  PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 \
  OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/bizon/anaconda3/bin/python -B \
  benchmarks/work/publication_qfo_parameter_swiss_audit_v1/benchmark_tools/audit_qfo_parameter_swiss.py \
  --inventory benchmark_tools/results/qfo_parameter_control_count_inventory_20260919.json \
  --inventory-sha256 25fbd435e7afb1dddd410433dc6df343c3c4586f3fa81a45980c39a19cdefa8c \
  --baseline benchmark_tools/results/qfo_swiss_counts_20260917.json \
  --output benchmarks/work/qfo_parameter_control_counts_20260919.json
```

Log: `benchmarks/work/qfo_parameter_control_audit_21953.log`.
All six changed variants remain explicitly not admitted. The check cannot
produce variant contrasts or stand in for the outstanding CPM workflow.
Terminal accounting and output validation must precede retaining its result.

## Completed Result

Job 21953 completed successfully, exit 0:0, in 12 seconds with 2 CPUs on
bizon. The batch log is empty and the result JSON is complete. Retained
`qfo_parameter_control_counts_21953.json` (245,132 bytes), SHA-256
`89ead915b35c90a00daa9664bbef6e0a01db2cf9c515535fd336b4ea1ec2cac7`.
The retained copy's hash equals the scheduled output.

The audit checked 733 file identities and reconstructed all 18 families
over 10,765 reference relations. Control F1 is 0.8335132180095781,
precision 0.9551767493219016 and recall 0.7393412564640668. A separate
structured comparison against the completed corrected factorial count
artifact confirms exact equality of the entire control family inventory,
raw counts, represented genes, family statistics and aggregate, not merely
rounded F1. The original reference relation count also matches.

All six variants remain `not_admitted`. Both `uncertainty_admitted` and
`publication_ready` remain false. This verifies the new auditor's real-data
control integration; it does not complete parameter sensitivity analysis or
produce confidence intervals.
