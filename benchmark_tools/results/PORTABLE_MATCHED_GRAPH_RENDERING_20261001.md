# Explicit Portable Matched-Graph Rendering

Previous goal turn is progress at pushed `9b3bd5fb`: test-only security pins,
511 focused checks and actual one-ULP macOS diagnostic evidence. Reread the
full objective and ledger. CI 36949483608 is now terminal failure; no restart.
At the recorded 01:34:36 UTC snapshot, source-9b3bd5fb run 36950764034 has
successful Linux CPU-wheel/docs jobs and five still-running test jobs.
Do not infer their outcomes or repeat the earlier wheel fixture.
At the later 01:40:22 UTC observation this existing run is terminal failure:
all five test jobs failed, Linux wheel/docs succeeded. No new sibling log is
inspected here, so their individual causes are not inferred. Do not restart it.

## Numerical Contract

The [actual inspected macOS difference](ci_exact_replay_difference_20261001.json)
is one ULP in an exploratory recall interval, not a changed plotted F1 or
bootstrap setting. Its kernel cause remains unisolated. The
[26 September standalone count protocol](MATCHED_GRAPH_STATISTICAL_REPRODUCTION_20260926.md)
already independently reproduces all 24 effects, eight adjusted F1 intervals
and seed signs within absolute 1e-12, relative zero. That contract predates
this rendering failure; no tolerance is selected by searching for green CI.

Keep renderer default `exact`, including one-ULP rejection and exact engine
metadata checks. Add explicit opt-in `count-level` through both API and CLI.
There is no fallback on exact failure and no user-selected tolerance. Portable
mode requires the unchanged independent count reproducer to pass, then compares
the complete contrast/bootstrap structure. Only finite floating contrast leaves
within that existing absolute tolerance may differ; metadata, structure,
missing/extra fields, inference roles and differing discrete counts still reject.
Every accepted exact difference and the complete count-validation report are
written to the export manifest with the policy, tolerances and scorer source.

The chart uses the retained reported values, not replacement statistics.
Historical results, figures, source bundles and review PDFs are not overwritten.
Exact identity and portable numerical reproduction are different claims.
NumPy metadata remains checked; this does not allow arbitrary engine versions.

```sh
python -m benchmark_tools.plot_matched_graph \
  --results benchmark_tools/results/matched_graph_scores_20260926/results.json \
  --output /absolute/fresh/matched_graph_portable \
  --replay-policy count-level
```

Leaving off `--replay-policy` retains strict behavior. Output must be fresh.
Record the input before reading and check it again after export; changed input
cannot receive a success manifest. Retain failed outputs rather than hiding
them, close figures on export failures, and refuse an existing destination.
Manifests continue to state `visual_review_complete: false`; an export is not
a new visual review or native/scientific workflow admission.

## Executed Checks

In the existing separate updated-test environment, **101 focused cases pass
in 7.72 seconds**, zero failures/errors/skips. Forty-nine new cases are included:
explicit exact-current-runtime acceptance, preserved one-ULP rejection, the
observed remote difference as a synthetic admission fixture with actual count
reproduction, all 24 effect and eight adjusted-interval corruption checks,
metadata/structure/nonfinite/count/panel rejection, independent-count enforcement,
unknown-policy rejection and actual CLI export/input-stability checks.

The retained-report chart test now explicitly tests portable rendering. A
separate strict test uses freshly calculated current-runtime values; it does
not silently modify a historical report or reinterpret exact cross-platform
failure as success. The CLI test really exports PNG/PDF/SVG into its temporary
fixture and verifies every output/input record and manifest policy. No native
inference or original data is opened by these count computations. Counts overlap
prior panels and are not additive; the earlier 93-case snapshot is superseded.
The [receipt](portable_matched_graph_rendering_20261001.json) pins JUnit and
all relevant code/evidence; scorer, standalone reproducer, retained result and
application/test requirements remain byte-equal to the preceding commit.

## Security And Remaining Work

A separate post-push read-only API check confirms
[all fifteen current test-manifest alerts are `fixed`](test_dependency_alert_closure_20261001.json),
none dismissed. Forty open alerts remain on four historical locks; those
locks are unchanged and are not recommended fresh installations. This is
repository alert-state evidence, not comprehensive runtime/native/OS security.

Commit/push the tested explicit policy, then inspect its actual automatic CI
handle without resubmitting earlier runs. New macOS execution and broader
failures remain open. Timing stays deferred without host polls/questions,
DGX, shared-package changes or unrelated job/service work. Uncertainty/source
gaps, runtime/data rights, controlled resources and the public release remain
incomplete; the full publication goal remains active.
