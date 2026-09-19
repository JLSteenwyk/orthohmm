# Parameter Uncertainty Integration

## Implementation

Frozen executor `1975265c4f58d24afe52617280d4843c8be53905` is pushed and
retained at `benchmarks/work/publication_qfo_parameter_uncertainty_v1`.
This extends the existing corrected SwissTrees count auditor to the two
prespecified CPM arms, then binds reconstructed counts to the existing
six-contrast numerical bootstrap. The original control and norm/margin
routes remain intact. Already queued executors are not modified.

CPM count admission requires the separately frozen score-admitter hash,
matching arm/index/context, native-pair conversion with zero mapping loss,
participant identity and matching execution evidence. Raw relation labels,
members, orientation and all family confusion counts are reconstructed and
checked against admitted native scores. Historical prediction counts are
not reused; the pinned historical audit supplies only reference truth.

The new runner pins the count auditor, numerical kernel, relevant validators,
protocol and plan. It always uses 100,000 shared PCG64 family draws, seed
20260925, linear percentile intervals and the full 18-endpoint multiplicity
denominator. Unavailable arms remain explicit and have null metrics. A
missing control makes every contrast unestimable. Partial panels cannot be
reported as complete, and no contrast is admitted from a control-only run.

Source SHA-256:

- `audit_qfo_parameter_swiss.py`:
  `f0b8d9256febb97041a6c1dd960cc42837ada14a69e1d1d03a921b7ab059eef3`.
- `run_qfo_parameter_uncertainty.py`:
  `6d1926ea886daefd5c8150630d9e4a03aa904daed4f0e01b0de2821d0df5ed44`.
- Control-only inventory:
  `2cad3c5d168c49b4337668f14559f001dc1b7d7cd18b83f69668a208260e88c9`.

## Real-Data Integration Check

Slurm job 21986 completed exit 0:0 in 18 seconds on bizon, using 2 CPUs
and a requested 32 GiB, one-hour limit, no requeue. Terminal accounting was
checked before inspecting output. The batch log is empty. This is an
integration check, not a new scientific parameter result or timing benchmark.

```bash
sbatch --parsable --job-name=qfo_parameter_uncertainty_control \
  --nodelist=bizon --cpus-per-task=2 --mem=32G --time=01:00:00 --no-requeue \
  --output=benchmarks/work/qfo_parameter_uncertainty_control_%j.log \
  --wrap='env -u PYTHONPATH -u PYTHONHOME -u LD_PRELOAD -u LD_LIBRARY_PATH PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/bizon/anaconda3/bin/python -B benchmarks/work/publication_qfo_parameter_uncertainty_v1/benchmark_tools/run_qfo_parameter_uncertainty.py --inventory benchmark_tools/results/qfo_parameter_uncertainty_control_inventory_20260919.json --inventory-sha256 2cad3c5d168c49b4337668f14559f001dc1b7d7cd18b83f69668a208260e88c9 --baseline benchmark_tools/results/qfo_swiss_counts_20260917.json --plan benchmark_tools/results/qfo_parameter_neighborhood_plan_20260919.json --protocol benchmark_tools/results/QFO_PARAMETER_NEIGHBORHOOD_PROTOCOL_20260919.md --output benchmarks/work/qfo_parameter_uncertainty_control_20260919.json'
```

Retained result: [qfo_parameter_uncertainty_control_21986.json](qfo_parameter_uncertainty_control_21986.json),
462,405 bytes, SHA-256
`c930b10fde59c799bc8d7f297a958320d07e5b45c934433b6afe49cafc837c47`.
The copy matches the scheduled output hash. The report lists 754 checked
file records, 18 families and 10,765 reference relations. A structured
comparison confirms exact equality of the reconstructed control arm,
including raw-file identity, all family counts, represented genes,
per-family statistics and aggregate, against the prior control audit 21953.

| Control endpoint | Value |
| --- | ---: |
| F1 | 0.8335132180095781 |
| Precision | 0.9551767493219016 |
| Recall | 0.7393412564640668 |

All six comparisons are `not_estimable`, with null metrics and family
differences. `estimated_contrasts` is zero; `uncertainty_admitted`,
`complete_panel` and `publication_ready` are false. The inventory explicitly
states that CPM admission 21984 and norm/margin admission 21948 are queued.
It is not subsequently mutated to incorporate later outcomes.

## Verification and Next Actions

199 focused tests passed in 3.62 seconds across the updated auditor,
provenance wrapper, real numerical kernel, CPM scoring/admission and shared
family-statistic helpers. Cases cover both CPM arms, wrong context/source,
boolean indexes, changed truth/member inventories, missing arms/control,
fixed multiplicity, changed protocol/input/helper hashes, mid-run mutation,
and output overwrite refusal. Frozen CLI import and staged whitespace
checks passed. This is not a fresh full-repository test run.

After genuine score admissions become available, freeze a new complete
seven-arm inventory, retain missing/failed/pending outcomes explicitly, and
run this executor on those records. Independently reproduce numerical
contrasts and generate the endpoint table/figure. Do not promote a default,
reduce the multiplicity denominator, or treat a confidence interval that
contains zero as evidence of equivalence. Independent generalization,
matched-search sensitivity and controlled scaling remain separate gaps.
