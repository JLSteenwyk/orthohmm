# First QfO Factorial Reconciliation Validated

The profile-refinement-off, candidate-expansion-off, reconciliation-on cell
`p0_c0_r1` completed as Slurm task21671_0 (raw job21672): COMPLETED0:0,
32CPUs on bizon, scheduler elapsed1:14:06. This elapsed time includes the
execution wrapper and is not a dedicated-machine scaling measurement.
Initial HMM search remains present in this cell.

## Independent Validation

Executed the existing frozen admission script from executor
`adec7e1010d1ab87424064db181600d0e44b7dfc` before the full reconciliation
array finished. The separate early-review output does not replace the
scheduled admission reports or bypass conversion/scoring gates:

```sh
python benchmarks/work/publication_qfo_factorial_admission_v1/benchmark_tools/admit_qfo_factorial_cell.py \
  --root "$PWD" --index 0 \
  --output benchmarks/results/qfo_factorial_v1/native_admission_0_early_review_20260917.json
```

Committed snapshot: `qfo_factorial_native_p0_c0_r1_20260917.json`.
The verifier checked successful scheduler identity, exact frozen execution
and input provenance, runtime tools, 157,073 recorded output artifacts,
native completion records, species-tree coverage and finite branch lengths,
and complete partition coverage. All976,504 candidate genes are retained
in399,387 root HOGs from393,231 candidate families;1,568 source families
split and no root HOG merges distinct source families.

All4,966,346 native cross-species ortholog pairs are canonical, unique,
sorted, within their source candidate family, and consistent with the native
summary count. Their412,214,095-byte file has SHA-256
`fdb0d07cd0924405b44f3aa326850f5050615232b4652ca1e0ddfca6aa5fd90f`.
QfO must evaluate these native pairs, not root-HOG clique pairs.

The native manifest reports26,171 reconciled and367,060 bypassed families.
Its internally inferred species tree uses17 families and20,459 supermatrix
columns, with one placement family. These are provenance observations, not
an independent test of species-tree or duplication-event correctness.

## Boundaries And Next Steps

This validates native output integrity, not orthology accuracy. Reference
mapping, conversion validation and native QfO assessment remain outstanding.
The scheduled admission/conversion arrays remain unchanged and will repeat
their own gates. Three of eight factorial cells currently have admitted
scores, all reconciliation-off; no missing score is imputed.

The existing admission unit suite passed17tests. No inference settings,
endpoints or statistical protocol changed, and no scientific run restarted.
The complete eight-cell comparison and42-endpoint SwissTrees analysis remain
pending, as do the broader publication requirements.
