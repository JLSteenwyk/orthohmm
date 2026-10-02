# Frozen Simulation Arithmetic Replay

This standalone workflow reproduces arithmetic from the two admitted native
simulation reports. It imports no project modules and reads no native outputs,
truth memberships, simulator inputs or historical paths retained in the JSON.
It does not repeat native inference, scoring or admission. Settings and results
stay frozen; this is not a complete executable study release.

## Inputs and Environment

Use Python 3.12.3 and NumPy 2.2.6 for the recorded replay. NumPy is the only
nonstandard import; its compiled/OS runtime is not bundled or asserted hermetic.
Use a private environment rather than upgrading a shared scientific interpreter.
The existing `requirements-ygob-arithmetic.txt` pins the same numerical dependency.

| Input | Bytes | SHA256 |
| --- | ---: | --- |
| `simulation_fixed_native_results_20260916.json` | 795,597 | `305dcb1dde0c0f57d8b390b0f00cc6148d95103a7efb98bd9f6dce37b96e64a8` |
| `simulation_variable_native_results_20260916.json` | 2,217,675 | `cc99fc31c3433809098212d3c8dc12f4ad829b4cfeb66dabcb13aede4e95258f` |

Transfer these already-retained reports and `reproduce_simulation_panels.py`.
Verify the script's identity independently before execution. No raw data,
native tools, Git checkout or new derived observation is needed for this replay.
Historical paths in the reports remain provenance only, not required inputs.

```bash
python -I -B /copied/reproduce_simulation_panels.py \
  --results /copied/simulation_fixed_native_results_20260916.json \
  --results /copied/simulation_variable_native_results_20260916.json \
  --output /fresh/simulation-arithmetic.json
```

Choose an unused output path with an existing parent. The script requires both
distinct panels, exact input bytes/digests and NumPy version, bounds report
size, rejects symlinked inputs, checks input bytes before/after replay and
refuses output overwrite. Validation failure produces a separate failed receipt.

## Statistic and Failure Boundary

Per-seed F1, precision and recall are independently computed using integer
TP/FP/FN and rational arithmetic. Pair-universe identities, endpoint coverage,
undefined ratios, truth sharing and diagnostic checkpoint/full-parent admission
are checked. All 280 planned outcomes per panel and their reasons/exclusions
remain explicit. Failures are never assigned zero accuracy.

Method means average available seed metrics; paired contrasts use only complete
paired seeds. Gene pairs and the two length panels are not pooled. Bootstrap
replay uses the original 20,000 PCG64 multinomial draws, reset per contrast,
with panel-specific seeds 20261031/20261130 and linear quantiles. Explicit
per-seed accumulation is independent of the producer's matrix product but uses
the same NumPy RNG/quantiles, not a new statistical engine. Values must match
within absolute tolerance 1e-12, without array broadcasting.

| Panel | Complete | Failed | Estimated contrasts | Unavailable contrasts |
| --- | ---: | ---: | ---: | ---: |
| Fixed length | 134 | 146 | 0 | 14 |
| Variable length | 267 | 13 | 14 | 0 |

Fixed-length OrthoFinder runs are not admitted: no paired effect or interval
can be inferred from that panel. Variable-length replay checks 42 paired metric
effects and 112 scalar interval bounds: nominal intervals for all three metrics
and 14-contrast Bonferroni intervals for primary F1 only. Precision/recall stay
exploratory; the sequence-only checkpoint is diagnostic, not an added contrast.
All 168 method-mean cells, including unavailable `null` values, are checked.

Native admission, generator correctness, truth history, exchangeable simulation
seeds and failure-related complete-case bias are not independently established
by replaying counts. Ten planned seeds limit interval resolution. This workflow
does not validate tree-perturbation or matched-recall panels, controlled timing,
arbitrary-dataset generalization, raw-data rights or full release readiness.

The [execution receipt](results/SIMULATION_ARITHMETIC_REPLAY_RESULT_20261002.md)
records actual copied-tree reproduction, failed harness attempts and tests.
