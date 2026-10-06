# Native SwissTrees Exclusions Localized

## Actual Finding

The admitted R1 input-cluster checksum equals the admitted R0 clustered-file
checksum. Reconstructing every source family from the admitted R1 RootHOG
records independently reproduces those exact candidate input bytes:
`4dadc740f64cd0a539a099e0577094e47ad4ba5862374dd8a8aca86a78eabdbf`.
The reconstruction covers394,328 source families,397,041 final RootHOGs and
all984,137 input genes. This demonstrates matching candidate partitions,
not identical upstream search-hit histories.

The [complete localization ledger](native_qfo_swiss_reconciliation_pairs_20261006_v1.tsv)
traces all2,023 changed reference relations from the preceding
[pair-decision analysis](NATIVE_QFO_SWISS_PAIR_TRANSITION_RESULT_20261006.md).
Every pair is present in a shared input candidate family, absent from the
admitted final native prediction file, and excluded at a positive-paralogy
duplication node with child-species overlap. The check reads all5,113,820
native prediction rows,858,427 annotation rows, and validates topology and
event-rule consistency for2,248 nodes in all23 affected source families.

All1,689 removed false positives remain within the same final RootHOG.
Of334 removed true positives,280 remain within a RootHOG and54 cross final
RootHOGs. The latter comprise28 RPS and26 VATB relations. All211 removed
PSEN true positives remain within a RootHOG. Consequently, on this observed
changed-pair universe, final group splitting alone does not account for the
precision gain or most recall loss: within-group pair resolution excludes
the relations even when group membership is unchanged.

The observed annotations label1,688 removed FP and332 removed TP exclusions
as high confidence, and one FP/two TP exclusions as medium confidence.
These are the method's serialized model labels, not calibrated orthology
probabilities or evidence that the underlying inferred duplications are true.
Do not use these observations to tune the frozen method or claim calibration.

## Independent Check

The [saved-Newick readback](native_qfo_swiss_reconciliation_newick_readback_20261006.json)
uses Biopython1.87 rather than importing the primary node parser, parent-walk
or localization functions. It checks all23 saved family checkpoints against
raw/rooted/annotated-tree digests, identical leaf universes and identical
rooted/annotated topology. Its common-ancestor API independently matches
each pair's exclusion-node ID, duplication event code and descendant count.
It checks1,139 leaves,51 distinct exclusion LCAs and every2,023 localized
pair, retaining the exact original changed-pair inventory and summary.

Original RootHOG, manifest and native prediction records were inventoried by
the original scientific output validation. The global node TSV and selected
checkpoints/trees were not. They are newly observed current-byte evidence,
explicitly not retroactively admitted artifacts or execution traces. Hash and
topology consistency strengthen the software-path explanation, not the
biological correctness of inference. All direct supplied/selected evidence is
checked before and after its respective analysis; no tree inference, rooting,
reconciliation or endpoint scoring is rerun.

The primary [result](native_qfo_swiss_reconciliation_trace_20261006_v1.json)
has SHA256 `53c692cf7dab9affac7b0240cf5191d2530bd4369b66f4165e8acb811391dacb`.
The complete localization ledger has SHA256
`4ac42d9537dd6b914aa019b2b45005b393cb96d56ebdd941981509f4430d60a9`.
The independent Newick readback has SHA256
`ee83382000bb48b95323209a197e7edb44de918bad6c6b8a1fbf8a1ea6e38835`.

281 joined tests pass in12.54s, no skips, including28 primary topology/
identity/event-rule cases and23 actual Newick/checkpoint/corruption cases.
Coverage includes prior pair-transition/SQL/native count/uncertainty/export/
functional diagnostics. The real trace uses original Python3.10 and finishes
in15.55s with66,048KiB maximum RSS; these are shared-host postprocessing
observations, not native inference timing. Contention remains unknown and
potentially tool-dependent. Failed native R1 timing stays ineligible.

## Reproduction And Limits

The [retrospective protocol](NATIVE_QFO_SWISS_RECONCILIATION_TRACE_PROTOCOL_20261006.md)
and scripts retain all exclusions without selecting successful families.
From the repository root, use fresh output paths:

```bash
env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  benchmark_tools/trace_native_qfo_swiss_reconciliation.py \
  --transitions benchmark_tools/results/native_qfo_swiss_pair_transitions_20261006_v2.json \
  --transitions-sha256 fe003f2cbc4285ea56cd80a92b71703b5c9e135ad1c8f244504dcfeab0911b46 \
  --output benchmark_tools/results/native_qfo_swiss_reconciliation_replay.json \
  --ledger benchmark_tools/results/native_qfo_swiss_reconciliation_pairs_replay.tsv

benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  benchmark_tools/readback_native_qfo_swiss_reconciliation.py \
  --report benchmark_tools/results/native_qfo_swiss_reconciliation_trace_20261006_v1.json \
  --report-sha256 53c692cf7dab9affac7b0240cf5191d2530bd4369b66f4165e8acb811391dacb \
  --output benchmark_tools/results/native_qfo_swiss_reconciliation_newick_replay.json
```

The second command independently checks the retained result, not the newly
created replay. Retained raw/native artifacts must remain locally available;
this is not whole-study archival restoration. Existing manuscript/PDF/archive
bytes remain unchanged; incorporate the companion evidence at the next
complete manuscript assembly.

This establishes a measured software exclusion path, not true duplication
history, inferred-tree accuracy, independent generalization, all-prediction
coverage or selected-default superiority. Initial HMM search is on in both
cells. Direct search evidence, prespecified error strata and tree-error
analyses remain separate work. No new confidence intervals or scientific
admissions, next-native authorization, inference retry or timing repair.
The native index8 workflow remains independent and pending; the full
publication goal remains active and unproven.
