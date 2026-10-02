# Current Accuracy Overview Reconciliation

The [new overview PDF](figures_current_accuracy_overview_20261002_v2/current_accuracy_overview.pdf)
and [24 plotted values](figures_current_accuracy_overview_20261002_v2/plotted_values.tsv)
now agree with the unchanged
[current eight-method score table](current_benchmark_scores_20260926_v2/scores.md).
The older overview used historical SonicParanoid Three Kingdoms F1 (0.9907944085)
instead of its later matched-input result (0.9912758997). That is a presentation
difference of +0.048149 percentage points, not a new method improvement. All
other Three Kingdoms points agree within floating-point roundoff. Historical
figures, predictions, scores, native admissions and defaults remain unchanged.

## Verification And Interpretation

The new exporter reads six hash-pinned direct records: the current score table,
corrected QfO comparison, matched-input Three Kingdoms comparison, OrthoBench
provenance register, historical overview manifest and corrected QfO figure
manifest. It verifies both frozen plotting helpers and checks all 72 score-table
cells, including the project-defined secondary mean, against their source
records. Three Kingdoms precision/recall/F1 are also recomputed from saved
TP/FP/FN counts. This is sufficient-statistic arithmetic, not raw biological
rescoring, conversion or new native admission.

All 48 coordinates and output semantics in the
[existing corrected QfO figure](figures_corrected_qfo_endpoints_20260926/corrected_qfo_endpoints.pdf)
match its complete source, and all three PNG/PDF/SVG artifact hashes still
match. That figure is not rerendered. Its PNG was inspected again: axes, legend
and notes are readable, while near-equal points can overlap. The exact table
and manifest, not visually separated glyph counts, identify all 48 points.
Its PDF was not newly visually inspected in this checkpoint.

The overview preserves all eight methods in fixed order. OrthoBench displays
weighted group precision/recall; Three Kingdoms displays projected BUSCO-reference
pair F1, including within-species pairs. Non-reference-gene errors are outside
that latter endpoint. Only SonicParanoid has the contemporary matched-input
row; other historical consumption gaps remain. This is not a newly matched
eight-tool experiment. OrthoFinder's sequence row is a pre-phylogeny MCL
checkpoint, and FastOMA uses supplied trees. GO/EC/FAS are similarities, not
F1. No cross-dataset average, universal ranking, new confidence interval or
controlled resource comparison is defined.

The extended manuscript now labels its original-input QfO comparison explicitly
as historical and links this current overview separately. The condensed main
text remains unchanged and already uses corrected QfO values. Its latest
eight-page render remains a dated snapshot, not a rendered copy of these later
extended-source edits. No new manuscript render or archive is generated here.

## Visual Review And Reproduction

The first new overview export failed visual review because an added caption
overlapped the Three Kingdoms x-axis label. Preserve that attempt locally and
move only the caption into the lower margin. The first/final plotted TSVs are
byte-identical. Final PNG and a bitmap of the one-page PDF were inspected;
labels, notes and legend are readable with no observed clipping or overlap.
PDF block-bound checks have zero violations. The
[review receipt](current_accuracy_overview_review_20261002.json) retains all
input/output pins, the earlier failed-visual artifact hashes and test evidence.
The final PDF is 23,460 bytes, SHA-256
`00d9e627c1c1a2ddb9f34164df30a71cebc2ff62c0be1651e88c6c719d4bfea8`.

```bash
python -m benchmark_tools.plot_current_accuracy_overview \
  --results benchmark_tools/results --output /absolute/new/current-overview
python -m pytest tests/unit/test_plot_current_accuracy_overview.py \
  tests/unit/test_plot_corrected_qfo_endpoints.py \
  tests/unit/test_export_current_benchmark_scores.py
```

The plot/export/source-table panel passes 50 cases in 3.27s before the final
prose check. It covers wrong scores/counts/units, inventories, provenance,
semantics, coordinates, means, source hashes and output preservation. The new
caption regression checks actual renderer boxes against axis labels/legend.
Passing tests did not replace the visual review that exposed the first issue.

The first extended-prose check failed on an exact limitation phrase: 59 pass,
one failure in 3.66s. Clarify the sentence without changing numbers or weakening
the assertion and preserve that failed receipt. The final combined panel,
including nine current-main-prose checks, passes **60 cases in 3.65s**, zero
errors/failures/skips. Its XML is 9,214 bytes, SHA-256
`ca7945fbb17deafc705b32a8e94167a6e75d87425fd2d2ebd4eefb72c0559140`.

## Prior CI And Open Requirements

Preceding source-5e899502 run 37006238922 at 12:32:56 UTC has successful
docs/Linux/wheel jobs, three failed fast macOS jobs and live 3.11/full jobs.
Its actual terminal Python 3.13 job 110835074599 log, downloaded once at
12:35:19 UTC, confirms source and 14,124 passes/four failures/119 skips/
30 warnings in 400.21s. All 12 null-figure, nine current-main-prose and ten
historical-artifact cases pass. This does not remotely confirm the new overview;
four unconfigured raw-export tests still lack original inputs/options. No
restart, inferred sibling counts or passing full-CI claim.
At 12:41:01 UTC all five macOS test jobs are terminal failures; docs/Linux/wheel
remain successful. That later state does not supply uninspected sibling logs
or counts and does not turn this new overview into remotely verified evidence.

This reconciliation closes a current-figure mismatch, not all figure/runtime/
rights or study-restoration gates. Timing remains deferred without questions,
host polls, DGX or unrelated process/service actions. Other-QfO uncertainty,
controlled resources, complete executable/versioned release, deposition and
final manuscript/figure/archive reconciliation remain open. Publication readiness
is not established and the full goal remains active.
