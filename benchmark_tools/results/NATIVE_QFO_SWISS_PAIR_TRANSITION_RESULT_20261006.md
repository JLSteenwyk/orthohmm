# Native SwissTrees Pair-Decision Trace

## Actual Result

For the already-admitted native P0/C0 cells, all10,765 assessed SwissTrees
reference relations have identical pair identities and truth labels in R0 and
R1. R1 makes no new positive predictions on this universe. It removes1,689
false positives and334 true positives. All2,023 changed decisions are retained
in the [complete pair ledger](native_qfo_swiss_changed_pairs_20261006_v2.tsv),
not a selected set of examples.

| Raw assessed label | R0 | R1 | Change |
| --- | ---: | ---: | ---: |
| TP | 3,123 | 2,789 | -334 |
| FP | 1,708 | 19 | -1,689 |
| FN | 1,193 | 1,527 | +334 |
| TN | 4,741 | 6,430 | +1,689 |

These pooled, unsmoothed counts are descriptive, not the native prior-adjusted
macro endpoint. Native SwissTrees F1 remains0.6891839606644314 for R0 and
0.7895743384374164 for R1. The existing conditional family intervals remain
separate [evidence](native_qfo_swiss_uncertainty_binding_22449_20261006.json);
this diagnostic adds no bootstrap draws or uncertainty admission.

The complete18-family transition tables identify where decisions change.
PSEN loses211 true positives and removes58 false positives; RPS loses28 true
positives without removing any false positives. APP removes148 false positives
without losing true positives. NOX removes552 false positives and loses one
true positive. Every family's results, including unchanged ASTER, are in the
[machine-readable trace](native_qfo_swiss_pair_transitions_20261006_v2.json).
These are descriptive locations for future biological/event tracing, not
cherry-picked validation or proof of incorrect trees.

## Validation

The original Python3.10 runtime replays direct scientific snapshot/admission
bindings and checks both original count-audit sources, raw execution/inventory
membership, exact family counts/members, original reference, native macro
arithmetic and unchanged truth before computing transitions. Both transition
marginals are checked. The primary result is28,028 bytes, SHA256
`fe003f2cbc4285ea56cd80a92b71703b5c9e135ad1c8f244504dcfeab0911b46`.
The complete changed-pair ledger is50,460 bytes, SHA256
`6276f2c7610accfd32332b1f90fa6d3191d901e3b1aed09340a78cd7b1c731ea`.

The actual [independent SQL readback](native_qfo_swiss_pair_transition_sql_readback_20261006.json)
runs stdlib-only under `-I -S -B`, without importing the primary diagnostic or
its statistics functions. It confirms21,530 raw rows,10,765 matched relations,
all288 family-transition cells,2,023 exact changed ledger rows and both native
macro statistics. Its9,195-byte result has SHA256
`09249d97c6ed3da6c840586552b4f5a72480cc46bbb944bc60124f5382164ffb`.
It verifies the same retained data, not a new biological replicate.

230 joined tests pass in7.64s, no skips, including33 primary trace cases and21
SQL/actual-artifact corruption checks. Joined coverage includes original
ordinary/recovered count audits, native scientific exporter, conditional
SwissTrees interval binding and functional-pair diagnostics. Primary integration
fixtures explicitly stub snapshot replay; actual execution does not.

The initial diagnostic exits1 before scientific outputs because it assumes an
ordinary snapshot contains `timing_eligible`. Preserve that attempt's resource
receipt and earlier tests. The corrected diagnostic retains timing fields only
where actually present; the recovered R1 timing remains false/null/ineligible.
The freshv2 execution completes in0.60s with98,304 KiB maximum RSS. These are
shared-host diagnostic/postprocessing observations, not native inference times
or isolated performance estimates. Contention effects remain unknown and
potentially tool-dependent. No original native job or scorer was rerun.

## Reproduction And Scope

The [protocol](NATIVE_QFO_SWISS_PAIR_TRANSITION_PROTOCOL_20261006.md) records
the retrospective scope and failure handling. From the repository root, use
fresh output paths for replay:

```bash
env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  benchmark_tools/trace_native_qfo_swiss_transitions.py \
  --snapshot benchmark_tools/results/native_qfo_scientific_scores_20261006_v1/report.json \
  --snapshot-sha256 6b2f735ea8a44f72715e11a7f015e4328576baa74c544889cbed0cdfe70ea07b \
  --counts-audit benchmark_tools/results/native_qfo_swiss_family_counts_20261005_v1.json \
    4a49225b50de31f153c1b701f866db7acdfdbefa3d6398339ef5bf8711b2b561 \
  --counts-audit benchmark_tools/results/recovered_native_qfo_swiss_counts_22449_20261006.json \
    af4276cd759412a5a10ba04012a6f742f147c28611a6cff83adee463f85b0471 \
  --output benchmark_tools/results/native_qfo_swiss_pair_transitions_replay.json \
  --changes benchmark_tools/results/native_qfo_swiss_changed_pairs_replay.tsv

benchmarks/work/native_factorial_review_py310_20261004/bin/python -I -S -B \
  benchmark_tools/readback_native_qfo_swiss_transitions.py \
  --report benchmark_tools/results/native_qfo_swiss_pair_transitions_20261006_v2.json \
  --report-sha256 fe003f2cbc4285ea56cd80a92b71703b5c9e135ad1c8f244504dcfeab0911b46 \
  --output benchmark_tools/results/native_qfo_swiss_pair_transition_sql_replay.json
```

Raw/submitted artifacts must remain locally available; this is not an archived
whole-study reproduction. The SQL command independently rechecks the retained
v2 result, while the first command makes a fresh primary replay.

This supports an exact pair-selection account of the observed trade-off, not
causal reconciliation attribution, tree-event accuracy, selected-default tool
superiority or independent generalization. Initial HMM search is on in both
cells and R changes clique to resolved-pair semantics. Further search/group/
tree-event tracing and other-endpoint uncertainty remain open. The pending
index8 workflow is separate; this trace does not admit its scores, authorize
index9 or establish publication readiness. Historical main-text/PDF/archive
bytes remain unchanged; incorporate this new companion evidence in the next
complete manuscript assembly after remaining native results are admitted.
