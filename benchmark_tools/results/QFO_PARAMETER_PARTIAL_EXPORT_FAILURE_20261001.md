# Retained Partial-Panel Export Failure

The first real-data export using `plot_qfo_parameter_neighborhood.py` at
commit `a06ae740` exited 1 before creating its destination. No figure, table,
new accuracy admission or uncertainty result was produced. Do not silently
retry or overwrite this attempt.

## Attempt

From the repository root:

```bash
/home/bizon/anaconda3/bin/python -B benchmark_tools/plot_qfo_parameter_neighborhood.py \
  --results benchmark_tools/results/qfo_parameter_uncertainty_partial_20260923.json \
  --sha256 1d6bd4f68a9d6728e0180bdde31f55ca0da7b4fc22e09bcc797dcb8db6cec628 \
  --reproduction benchmark_tools/results/qfo_parameter_uncertainty_partial_reproduction_20260923.json \
  --reproduction-sha256 252c8c0b16431c2b4993bb53b360d73405776555cbe288f70af9339881240774 \
  --output /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmark_tools/results/qfo_parameter_partial_export_20261001
```

The exception was `ValueError: Frozen input/source identity changed:` for the
root `benchmark_tools/export_qfo_corrected_comparison.py`. The traceback
reached the first `check(ref)` loop, before destination creation. Follow-up
filesystem checks confirmed that the destination did not exist.

## Diagnosis

The retained report has 5,125 checked-input records, representing 5,108 unique
paths with no conflicting identities. A targeted examination of its
`benchmark_tools/*.py` records found exactly one mismatching source:

| Source | Bytes | SHA256 |
| --- | ---: | --- |
| Retained root comparator exporter | 8652 | 4239f8d263295b313c5cc27b866bf7c6494feac44a8fbfe18e160bd92ff69bfb |
| Current root comparator exporter | 9572 | 8fd211aeaaac5765b6f5ffa01812d8e61ac4134b8c1ceb1bada4f497df7d7aa0 |

`git show 5efb206b23a44f85386d8cb7e90f3e815a3d162a:benchmark_tools/export_qfo_corrected_comparison.py`
produced exactly the retained bytes and SHA256. Comparing that revision with
`01ac4f66b2ac4b3565e5617b7abfe789e7ffe9aa` shows additions for recovered
OrthoMCL reporting, not parameter scoring or bootstrapping. The exporter was
transitively captured outside the explicit frozen scientific helper list.
This targeted examination is not full current-data admission.

The original exporter copy in the retained
`publication_qfo_cpm_candidates_admission_v3` worktree was rechecked: direct
path, 8,652 bytes, original SHA256. That worktree's HEAD is
`497c54f1de8010482c0aaabd92c2c4e68dff43da`; this observation alone does not
claim the whole worktree is clean. The current exporter also matches the
Git blob at `01ac4f66`. An explicit source-lineage route must retain both
historical and current identities and check every other binding unchanged.
Do not revert current sources, edit historical receipts or silently drop pins.

## Interrupted And Resumed Execution

Three continuations could not execute commands or save documentation because
the sandbox failed before process launch with
`bwrap: loopback: Failed RTM_NEWADDR: Operation not permitted`. The attempted
failure note and ledger edits did not succeed. No source-binding implementation
or successful export was claimed during that interruption.

On command-access restoration, assessment `22393` and its batch step were
COMPLETED, exit 0:0, 8 CPUs/64 GiB, bizon, elapsed 42:55. Its execution report
is `process_succeeded_pending_independent_admission`. Independent score gate
`22394` was submitted once from the existing clean frozen executor; its actual
submission is recorded separately. No current score is admitted by this note.
Timing remains deferred; no DGX or unrelated job/service action.
