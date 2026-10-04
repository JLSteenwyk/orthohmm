# Prepared Resource Reporting Component

The [new component](bundle_shared_prepared_component_20261004.py) archives the
prepared-panel table/figure format, including measured exclusions and
pre-native aborts with absent endpoints. It reuses the actual reporting and
plotting functions rather than independently implementing their arithmetic.
The historical v1 component/source/archives remain unchanged. New bundles
explicitly use schema 2 and include the hash-pinned historical standard-library
helper alongside the current prepared table/figure sources.

## Build And Replay

Use committed source bytes and a completed snapshot/figure pair. The component
can replay partial reports, but that does not make the 27-attempt panel complete.
For the final release, select the final reviewed pair, not this interim example:

```bash
python -B benchmark_tools/results/bundle_shared_prepared_component_20261004.py build \
  --repo . \
  --snapshot benchmark_tools/results/threadripper_shared_panel_snapshot_20261004_v22 \
  --figure benchmark_tools/results/threadripper_shared_resource_figure_20261004_v21 \
  --output /absolute/fresh/prepared-resource-component
```

Retain the emitted `bundle.json` SHA-256 outside the bundle before transfer.
After extraction, verification uses only standard-library Python and the
included readers. It checks the externally anchored inventory before
recomputing report consistency; original evidence paths are metadata only.

```bash
python3 -I -S -B /relocated/component/component.py verify \
  --directory /relocated/component --manifest-sha256 RETAINED_SHA256
```

For numerical/visual replay use the recorded Matplotlib version from
`requirements.txt` (3.10.8 for the current figure), with NumPy 2.2.6 and
Python 3.12.3 as used in the retained reporting environment:

```bash
/reporting-env/bin/python -I -B /relocated/component/component.py replay \
  --directory /relocated/component --manifest-sha256 RETAINED_SHA256 \
  --output /absolute/fresh/replayed-report
```

Replay compares all three table files byte-for-byte and rendered PNG pixels
exactly. PDF/SVG originals are checked and delivered, not claimed byte-identical
after rerendering. Use a fresh output outside the immutable component.

## Validation Scope

Eighteen focused tests pass in 3.97s. Tests use the actual 22-attempt table and
20-point figure: 19 eligible measurements, exclusions 0/17/20, pre-native aborts
17/20 and two complete cells. A fresh archive extracts and verifies through
the copied CLI with `-I -S -B`; exact tables/PNG replay. Negative cases reject
wrong anchors, changed/missing/extra/symlink members, scope/schema/coverage
promotion, imputed abort endpoints, altered medians and unsafe output reuse.
The unit builder mocks Git reads; separately building from the committed
repository is required for actual archive provenance.

This closes a reporting-format compatibility gap. It does not replay raw
accounting, scientific inputs, native inference, accuracy or uncertainty;
establish independent validation; remove contention; or complete the full
versioned publication release. Missing repeats and failures stay visible.
