# Controlled Fragment Workflow

This supplements the publication reproduction instructions for the fragment
observation and retained-stage diagnostic. It is not a new release candidate,
an all-tool runtime certification or evidence of publication readiness.

## Frozen Evidence

Run from the repository root. Keep these original results unchanged:

- `results/controlled_fragment_results_20261009_v1/report.json`: 80 baseline/
  fragment score records, 240 endpoint strata and 15 paired comparisons.
- `results/controlled_fragment_execution_20261009_v1.json`: native execution
  and independent score/interval readback commands and outcomes.
- `results/controlled_fragment_trace_selection_20261009_v1/selection.json`:
  84 bins, 53 selected unique records and 31 empty bins.
- `results/controlled_fragment_stage_trace_20261010_v1/report.json` and
  `stages.tsv`: 106 observations for those cases in both arms.
- `results/controlled_fragment_stage_readback_20261010_v2.json`: successful
  independent native-file/NA-table readback, not a new score or inference run.
- `results/controlled_fragment_trace_execution_20261010_v1.json`: actual
  selector, native tracer, independent readers and preserved failed readback.
- `results/controlled_fragment_trace_integration_20261010_v1/manifest.json`:
  exact input/source/output pins, 11 location-table rows and one explicit
  status-sentence replacement in the new manuscript.
- `results/controlled_fragment_trace_integration_execution_20261010_v1.json`:
  actual generation, independent table/parent readback, HTML/PDF execution,
  printed-cell check and scoped manual review.

Paths above are relative to `benchmark_tools/`. Original scientific inputs and
native output directories under `benchmarks/` are deliberately not committed
as large raw datasets. The prepared panel is
`benchmarks/work/controlled_fragment_observations_20261009_v1/manifest.json`.
The original full-length controls also use the prepared variable-length panel
and retained native executions. The selection's `bindings` identify every
baseline/fragment output, execution inventory and configured argument vector;
the stage report's `checked_inputs` identify exactly which files were read.
Those files must be retained to verify native stages. A clone containing only
Git summaries cannot reconstruct missing raw outputs or validate native causes.

Recorded paths are absolute to the original workspace. The frozen readers keep
those labels; `relocate_controlled_fragment_stages_v2.py` provides a separate,
read-only mapping to copied raw artifacts. Do not rewrite the old manifests or
claim that summary-only replay proves raw inference reproduction. The previous
relocated components retain their own documented, narrower scope.

## Relocated Native-Stage Check

The new component contains the stage report/table, selection and its pinned
inputs, expected independent readback, unchanged independent reader kernels,
and the small relocation runtime. Ownership FASTAs for every selected native
context are added only from retained execution input inventories; copied
OrthoFinder FASTAs for the sequence-only checkpoint use the parent's output
inventory. Preparation refuses changed inputs, conflicting mappings, an occupied
destination or more than 64 MiB of pinned payloads. Raw artifacts stay local;
this is not another release candidate or a public redistribution claim.

The first actual relocation attempt, source `43f1e398`, failed at ownership
lookup because 56 original configured FASTAs were absent from its copied map.
It and its component remain preserved in
`results/controlled_fragment_relocation_execution_20261010_v1.json`. The separate
`v2` planner reuses the frozen copy writer and native replay, adding ownership
bindings for all selected contexts. It records that preparation revision in the
same component layout, without changing native source or historical metadata.

The actual corrected component is retained locally at
`benchmarks/work/controlled_fragment_relocated_stages_20261010_v2`. Its 1,783
pinned files contain 39,428,185 payload bytes, including 112 copied-checkpoint
FASTAs and 56 additional configured originals. The
[retained map](results/controlled_fragment_relocated_manifest_20261010_v2.json)
has SHA-256 `cf42e0e5a9a0486879b976b5cc7848a94cb00d657caacbb64b77a1ba7d04f652`.
Python 3.12.3 with NumPy 2.2.6 executed the copied runner from a location outside
the checkout. The [actual readback](results/controlled_fragment_relocated_readback_20261010_v2.json)
matches the complete original result: 106 observations, 53 cases, 44 contexts
and 415 native checked-input references. The
[execution record](results/controlled_fragment_relocation_execution_20261010_v2.json)
retains preparation, copy, replay and independent object/source checks. A
Python-level `open` audit forbade reads of all 1,783 original artifact paths;
zero were attempted. This is a functional check, not OS containment. All raw
payloads remain local and are not supplied by a summary-only Git clone.

From the original repository, prepare at a fresh path:

```bash
benchmarks/work/release_alert_refresh_20261001/venv/bin/python -I -B \
  benchmark_tools/relocate_controlled_fragment_stages_v2.py prepare \
  --root "$PWD" \
  --report benchmark_tools/results/controlled_fragment_stage_trace_20261010_v1/report.json \
  --readback benchmark_tools/results/controlled_fragment_stage_readback_20261010_v2.json \
  --readback-sha256 2d37e635fb9d9f08d1b1857d829e0798c19452936785b92e5dd40a058658da35 \
  --output /fresh/path/fragment-component
```

Copy the complete component to its new location. Using the manifest SHA-256
printed by preparation, replay with Python and NumPy (no repository imports or
native inference executables are required):

```bash
python -I -B /new/path/fragment-component/runner/benchmark_tools/relocate_controlled_fragment_stages_v2.py replay \
  --component /new/path/fragment-component \
  --manifest-sha256 MANIFEST_SHA256_FROM_PREPARATION \
  --output /fresh/path/fragment-readback.json
```

The copied kernels are loaded privately and unchanged. Only their `Path`
binding is redirected; NumPy receives verified copied paths through
`os.PathLike`. Logical historical labels remain in results, with their physical
copies explicitly listed in `manifest.json`. Missing/unindexed originals are
never a fallback. Replay must reproduce the entire retained readback object,
including native numerical checks, summary and exact `NA` table serialization.
This verifies portable raw-stage postprocessing, not reproduction of historical
inference, benchmark scores, all-tool dependency closure or OS-level isolation.

## Execution Order

Source was committed before each production step: selector `9a9b076b`, retained
stage tracer `b2b2e222`, independent reader `e49fd970`, NA-only readback recovery
`49b808e5`, and manuscript integration `f10580cb`. Later documentation/results
commits do not change those attempt-bound sources.

1. The fixed fragment observation reuses seeds 20261101--20261110 and parent
   truth, truncating exactly the prespecified ranked subset. The prepared
   manifest and original execution receipt retain the exact commands/settings.
   All 30 native processes and 40 method/checkpoint outcomes completed. Do not
   automatically repeat them, the baseline controls or the timing panel.
2. `select_controlled_fragment_trace.run` selects minimum-hash records from
   all prespecified bins, including empty bins and retained TP controls. It
   executed once. Selection is retrospective, not independent generalization.
3. `trace_controlled_fragment_stages.run` reads already-retained stages. It
   executed once; no clustering, gene-tree inference or scoring is performed.
4. Use `readback_controlled_fragment_stage_table.verify` for the independent
   native readback. The first reader remains unchanged and failed only at TSV
   serialization: JSON null is exported as `NA`, not an empty cell or `False`.
   The corrected wrapper reuses independent native/graph/node kernels, not the
   stage producer. Never reinterpret `NA` as evidence that a hit was absent.
5. `integrate_controlled_fragment_trace.run` generates the new manuscript,
   section files, location table and bounded claims. Its manifest records all
   inputs and the explicit status replacement; removing the new sections and
   reversing that replacement restores the exact parent.
6. HTML uses `render_manuscript_review.render` and the existing scoped
   profile-table print correction. `print_manuscript_review.print_review` uses
   the normal browser sandbox. PDF checks and selected-page views are separate
   from numerical/native validation and do not certify all document pages.

## Dependencies And Checks

Actual native tracing and its independent check used the retained Python 3.12
scientific environment at
`benchmarks/work/release_alert_refresh_20261001/venv/bin/python`, with NumPy
for the independent checkpoint reader. The producer additionally uses BioPython,
SciPy and DendroPy-backed existing helpers. Native OrthoFinder MCL syntax parsing
is reused, but graph traversal and reconciliation/root reconstruction in the
independent reader are separate implementations. No HMM search engine runs in
these postprocessing commands. Original native scientific environments and
parameters belong to their execution receipts, not the current reader runtime.

Integration/HTML generation also ran with the existing minimal Python 3.10
environment. PDF printing/review used Python 3.12 with PyMuPDF 1.27.2.3; the
minimal interpreter's import refusal is retained and occurred before browser
launch. Chrome was `/opt/google/chrome/google-chrome`, without `--no-sandbox`.

Actual commands unset `PYTHONPATH`, `PYTHONHOME`, `PYTHONUSERBASE`, `LD_PRELOAD`,
`LD_LIBRARY_PATH` and `LD_AUDIT`; disable user-site imports/bytecode; and use
isolated Python. Numeric readers use one OpenBLAS/OMP/MKL thread. The exact
inline module-loading commands and terminal results are in the receipts above.
Production functions refuse occupied output paths. Reproduction must use a
distinct documentary output path and retain its source/input identity; it must
not overwrite prior attempts or silently trigger a new scientific run.

Focused tests, run using the retained scientific interpreter:

```bash
benchmarks/work/release_alert_refresh_20261001/venv/bin/python -I -B -m pytest \
  tests/unit/test_controlled_fragment_stage_readback.py \
  tests/unit/test_controlled_fragment_stage_table.py \
  tests/unit/test_controlled_fragment_trace_integration.py \
  tests/unit/test_fragment_trace_artifact_access.py \
  tests/unit/test_relocate_controlled_fragment_stages.py \
  tests/unit/test_relocate_controlled_fragment_stages_v2.py
```

These check native-reader kernels, unavailable-versus-false serialization,
bounded integration, exact parent restoration and changed-input refusals.
The relocation fixture also launches the copied standalone runner after moving
its component and removing all synthetic original files. It checks all four
methods, including numeric checkpoints and copied-input ownership. Tests alone
do not establish availability of the actual raw files on another machine.

## Interpretation

Counts describe selected representatives, not error prevalence or confidence
intervals. A connected graph is not orthology; a recorded duplication is not
proof of the inferred tree or biology. Missing significant hits do not isolate
prefilter versus scoring failure. OrthoFinder search/event adapters remain
unavailable, and its sequence-only checkpoint has no independent timing here.
Synthetic truncation does not validate natural fragment annotation or full
domain architecture. No initial-HMM-off benefit, parameter promotion, general
superiority or publication-readiness claim is made.

The completed timing panel remains shared-host evidence: contention has an
unknown, potentially tool-dependent impact. It does not require a quiet host
or DGX, and must not be repeated solely because other analyses were running.
