# Prepared Resource Reporting Component

The [new component](bundle_shared_prepared_component_20261004.py) archives the
prepared-panel table/figure format, including measured exclusions and
pre-native aborts with absent endpoints. It reuses the actual reporting and
plotting functions rather than independently implementing their arithmetic.
The historical v1 component/source/archives remain unchanged. New bundles
use schema 2 by default and include the hash-pinned historical standard-library
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
The unit builder mocks Git reads; the actual committed execution below
separately establishes archive provenance.

## Actual Committed Archive

Source commit `84217cc33fe111659940d852d8c30da75b926da3` builds the current
table v22/figure v21 once. The [execution receipt](prepared_resource_component_execution_20261004.json)
is 5,964 bytes, SHA-256
`4b2d1a84032704ecb3d3aee45f4999cf1fef9bc84284033bde347744e404934d`.
The local archive contains 19 payload files plus the manifest, all 20 regular
members checked for unique names/type/mode/bytes/hash before fresh extraction
under `/tmp`. Archive: 441,980 bytes, SHA-256
`7614c4b6c4ce5d48bbac13945ebe29e6ccaec0293033acf562fef096659920f3`.
External manifest anchor: 3,521 bytes, SHA-256
`edda5231a18e0634aa189dec6da03f8fd3c60468ee7d8e412a4462dbdfbc539e`.

The copied CLI verifies with `-I -S -B` and replays with `-I -B`, both zero
exit and no stderr. All three table files and PNG pixels match. Source/test/
historical-helper hashes remain unchanged; the retained JUnit has 18 passes,
zero failures/errors/skips. The archive remains under ignored local work paths,
not committed or uploaded. This is an interim reporting archive, not a final
panel or replacement for any older archive. No native audit or benchmark is
repeated to obtain it.

This closes a reporting-format compatibility gap. It does not replay raw
accounting, scientific inputs, native inference, accuracy or uncertainty;
establish independent validation; remove contention; or complete the full
versioned publication release. Missing repeats and failures stay visible.

## Prose-Inclusive Bundle Preparation

The optional `build --section /path/to/generated-resource-section` selects
schema 3 and adds the section generator, exact prose and its provenance
manifest. It checks the same table bytes, all five section-source pins,
coverage and reporting-only scope, then recomputes the prose. Relocated replay
now also emits byte-identical `resource_section.md`. This includes eight source
files and 23 regular archive members; it still does not open raw evidence paths
or reproduce native inference. Omitting `--section` preserves the schema-2
inventory. Historical archives and their bundled readers remain unchanged;
use their exact copied readers, not a newer reader substituted into old bundles.

Eight added tests plus the 18 existing component cases pass: 26 tests in 8.53s.
They use actual v24 data/figure v23/prose v24 with mocked Git reads and temporary
archives. Copied CLI verification runs with `-I -S -B`, copied replay with
`-I -B`; tables, PNG pixels and prose agree. Mismatched table/source/coverage,
scope promotion and reanchored edited prose are rejected. JUnit at
`benchmarks/work/threadripper_shared_execution_20261003/resource_section_bundle_tests_20261004.xml`
has SHA-256 `92e3d25acc0886c2545abce033bc429712a7e3ebbe6ff08d3d68788e29255da3`.
This is tested final-bundle preparation, not an actual committed-source final
archive execution. After all 27 attempt reviews, build the final component once
from committed sources and the final table/figure/section; independently
restore and inspect it. Do not rebuild unchanged interim archives per attempt.
