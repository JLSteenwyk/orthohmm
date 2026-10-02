# Frozen History And Exact Replay Diagnostics

The previous goal turn answered the user's timing question but made no new
execution or implementation progress. This continuation rereads the objective
and current ledger and verifies the actual existing CI handle, rather than
restarting it. Run [36947338909](https://github.com/JLSteenwyk/orthohmm/actions/runs/36947338909)
at `5b030d0fa1f9047da1d6ab89e0089233c5e01d83` is now terminal failure;
its independent Linux CPU-wheel and docs jobs succeeded, while all five macOS
test jobs failed. The preceding successful wheel report remains valid within
its recorded scope; it does not establish a passing full suite.

Download and inspect the actual Python 3.13 job log once. It reports
**13,364 passed, 201 failed, 47 errors and 96 skipped in 622.75 seconds**.
The runner is macOS 26 ARM64, Python 3.13.15; test installation installs
NumPy 2.2.6 after the application install. The inspected log still rejects
matched-graph rendering at the exact contrast/bootstrap guard but contains
no differing values. Do not infer floating-point drift, altered counts or a
metadata-only explanation from that message. The separate count-level
reproduction module passes in this log; it uses its already recorded 1e-12
tolerance, not the renderer's exact dictionary guard.

The [machine-readable receipt](ci_replay_diagnostics_20261001.json) pins the
retained ignored log, local JUnit and implementation/input identities. The
other four failing job logs have not been inspected here; their terminal
conclusions do not establish identical causes or test counts.

## Focused Changes

Two CI test checkouts now request `fetch-depth: 0`. The failing log contains
`git ls-tree`, `git show` and `git archive` errors when accessing scientific
revision `7f3a9e40dd7e79f842cc2c11fb8b548f9a802806`. These tests require that
actual historical source; substituting HEAD or skipping them would not test
the frozen archive. The [checkout documentation](https://github.com/actions/checkout#usage)
states that the default fetches one commit and depth zero supplies the full
history. Docs and CPU-wheel checkouts, test selection and all other workflow
settings remain unchanged. Remote confirmation is pending.

The matched-graph renderer retains its exact contrast/bootstrap comparison.
On rejection it now includes JSON paths and the reported/recomputed values,
including list indices and missing-key presence. This is failure diagnosis,
not an admission tolerance, changed statistics or a new bootstrap engine.
Existing plot inputs and scoring/count-reproduction code remain unchanged.
No historical rendered artifact or frozen result is overwritten. Remote
failure values must be inspected before deciding whether any numerical
comparison policy should change.

## Executed Checks

In the existing private Python 3.12 / NumPy 2.2.6 environment:

```sh
benchmarks/work/ci_numpy_pinned_20261001/venv/bin/python -m pytest \
  tests/unit/test_plot_matched_graph.py \
  tests/unit/test_score_matched_graph.py \
  tests/unit/test_reproduce_matched_graph_statistics.py \
  tests/unit/test_verify_frozen_source_archive.py \
  tests/unit/test_declared_test_dependencies.py \
  --junitxml=benchmarks/work/ci_numpy_pinned_20261001/matched_replay_guard_20261001.xml -q
```

**60 cases pass in 6.06 seconds**, no failures/errors/skips. Ten new diagnostic
cases are included: both one-ULP directions remain rejected, nested/list/shape/
missing-key diagnostics preserve exact values, and metadata mismatch is not
suppressed. The existing retained plot, count-level reproduction, corruption
rejection and actual frozen Git archive checks pass. Twenty-two historical
invalid-escape warnings from frozen source syntax checks remain visible.
These counts are not additive with prior overlapping test panels.

Structured YAML comparison proves only the two test checkout options changed;
all existing jobs, matrix settings, commands and unconditional CPU-wheel
artifact upload remain intact. Scoped whitespace passes. Local success uses
this checkout's full Git history and does not prove the remote fix has executed.

Commit and push this tested milestone, then observe the newly created CI
handle without manual resubmission. Accuracy, scientific defaults, historical
scores, controlled timing and published release status are unchanged. Timing
stays deferred without workload polling, quiet-window requests, DGX access,
shared-package changes or unrelated job/service actions. The complete
publication goal remains active and incomplete.
