# Retained Figure Integrity Update

The refreshed explicit inventory includes the completed **original-release**
QfO SwissTrees factorial figure. Its plotted results remain distinct from
the unfinished corrected-input factorial. The previous 14-panel audit is
retained unchanged as a historical snapshot.

The new [machine-readable audit](publication_figure_integrity_20260918.json)
verified all 15 manifests, 52 output records and 88 recorded file-reference
occurrences against saved byte counts and SHA-256 values. No changed or
missing recorded artifact was found. This was a local audit; no DGX files
were scanned during its active timing jobs.

The same single detached source dependency remains outside main-repository
tracking: `publication_method_native_v2/benchmark_tools/replay_high_sensitivity.py`.
Its frozen revision and byte recovery are documented in the
[previous audit](PUBLICATION_FIGURE_INTEGRITY_20260917.md). A final portable
bundle must export this source explicitly; its availability in a local
worktree is not a portability guarantee.

Ten focused audit tests pass, including explicit retention of the original
QfO factorial and exclusion of uncompleted corrected-factorial figures.
The separate [factorial statistical audit](QFO_FACTORIAL_SWISS_RESULTS_20260918.md)
and [figure manifest](qfo_factorial_swiss_figure_20260918/manifest.json)
remain the evidence for numerical provenance. The byte audit does not
independently establish plotted arithmetic, statistical validity, rendering,
transitive raw-data lineage, licenses or publication readiness.

```bash
python -m pytest -q tests/unit/test_audit_publication_figures.py
python benchmark_tools/audit_publication_figures.py --repo . --output /tmp/new-figure-integrity.json
```

Use a fresh output path. Matched scaling and corrected-input factorial
figures must not be promoted before their results pass the required gates.
