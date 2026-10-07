# Three-Cell Native QfO Figure

## New Presentation

The [new PNG](native_qfo_three_cell_figure_20261006_v1/native_qfo_three_cell.png),
[vector PDF](native_qfo_three_cell_figure_20261006_v1/native_qfo_three_cell.pdf)
and [SVG](native_qfo_three_cell_figure_20261006_v1/native_qfo_three_cell.svg)
present the three admitted fresh native cells P0C0R0, recovered P0C0R1 and
P0C1R0. Four score identities remain unavailable, not zeros or substituted
cached results. The existing two-cell figure/supplement and prior archive
payloads remain unchanged. This new presentation does not complete the
factorial, biological validation or publication package.

Panels A/B separate orthology F1 from GO/EC/FAS similarities. C shows
all-input relation coverage, not reference recall or accuracy. D/E show
candidate-expansion and reconciliation contrasts at P0, respectively;
both compare to P0C0R0. Thick intervals are nominal 95% paired-family
intervals; thin intervals retain the 42-endpoint adjustment. Both adjusted
F1 intervals include zero. F shows SwissTrees native precision-recall points.
No uncertainty is imputed for the other QfO endpoints; the project-defined
secondary mean is not plotted. Initial HMM search remains on in every cell.
R-off group-clique and R-on resolved-pair semantics remain distinguished.

## Checked Inputs And Outputs

[Manifest](native_qfo_three_cell_figure_20261006_v1/manifest.json) binds the
current scientific report `7916d3e...`, guarded three-cell SwissTrees binding
`76e6f4d3...` and independent candidate raw/rational readback `07941ba9...`.
The unchanged scientific binder replays the entire supplied binding exactly
in the original scientific Python3.10.13 environment (Bio1.87/NumPy2.2.6/
psutil7.2.2). Preserve its venv invocation instead of resolving its symlink.
Rendering uses the existing Python3.12 plotting environment, without changing
the frozen runtime, 14 allocation-route sources or 920 historical helpers.

Full-precision plotted data accompany the assets:
[18 scores](native_qfo_three_cell_figure_20261006_v1/scores.tsv),
[coverage and precision-recall](native_qfo_three_cell_figure_20261006_v1/coverage.tsv),
[six contrast endpoints](native_qfo_three_cell_figure_20261006_v1/swiss_intervals.tsv).
The [independent content readback](native_qfo_three_cell_figure_readback_20261006_v1.json)
checks those values exactly against the admitted snapshot/binding, output
hashes, SVG scope labels, actual one-page PDF text, dimensions and all three
cell colors. It decodes the actual PDF to a
[new preview](native_qfo_three_cell_pdf_preview_20261006_v1.png), rather than
accepting an unrelated preview file. Content/pixel checks do not certify
biological validity, interval exchangeability or visual layout automatically.

Separately inspect both the PNG and decoded PDF preview: all six panels,
legend, tick labels, coverage labels and scope footers are visible and readable,
with no observed incoherent overlap or clipping. The manifest retains
`visual_review_complete=false`; this document records the separate inspection
without rewriting the generated receipt. These checks do not prove all future
exports or viewing sizes are legible.

The initial [104-test suite](native_qfo_three_cell_figure_tests_20261006_v1.xml)
passes in 9.38s, zero errors/failures/skips. It covers preserved historical
assets, invalid scores, coverage denominators, null missing cells, failed-timing
status, both exact interval inventories, cross-input identities, no overwrite,
rendering and resealed table/label/blank-raster/PDF/source/hash failures.
An old documentation assertion is updated from historical two-cell counts
to the current three-cell claim; neither historical scientific source nor
figure receipt is changed to satisfy it.
The final [joined suite](native_qfo_three_cell_figure_joined_tests_20261006_v1.xml)
passes 362 tests in 15.07s, zero errors/failures/skips, including actual new
asset/source bindings and current manuscript/reproduction contracts. Actual
PNG is 2760x1960 (286,849 bytes); decoded PDF preview is 1656x1176. These
presentation tests are not a fresh audit of all raw biological evidence.

## Reproduce

Retained outputs already exist and are not regenerated on continuation.
Use fresh output, review and preview paths for a genuine reproduction.
Rendering is not inference, scoring, a new bootstrap or scientific admission.

```bash
env -u PYTHONPATH -u PYTHONHOME -u PYTHONUSERBASE \
  -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/release_alert_refresh_20261001/venv/bin/python -B \
  -m benchmark_tools.plot_native_qfo_three_cell_scores \
  --snapshot benchmark_tools/results/native_qfo_scientific_scores_20261006_v2/report.json \
  --snapshot-sha256 7916d3e23808edbb92b016b5c40ac56b9417a863d1af01a7b93a8ef5dbb53a63 \
  --swiss-binding benchmark_tools/results/native_qfo_candidate_swiss_uncertainty_20261006_v1.json \
  --swiss-binding-sha256 76e6f4d33a13a12b7f0a396f6b795b63b654666a038485fb8af2538259f553c2 \
  --candidate-readback benchmark_tools/results/native_qfo_candidate_swiss_readback_20261006_v1.json \
  --candidate-readback-sha256 07941ba9c55ba7dc14c5ba53df37d28cb3b2f60f599b9d177b10c9afa99c1329 \
  --validation-python benchmarks/work/native_factorial_review_py310_20261004/bin/python \
  --output benchmark_tools/results/native_qfo_three_cell_figure_20261006_v1
```

Review CLI: `benchmark_tools.review_native_qfo_three_cell_figure` with
`--figure-directory`, `--pdf-preview` and `--output`, in the same plotting
environment. Paths above identify the actual completed run, not instructions
to overwrite or duplicate its outputs. No new package is installed.

Next observe the same native identity-10 job23902; terminal review and
successful-output conversion/scoring/independent admission remain required.
11/12 stay sequential behind reviewed history; retain failure9 with no automatic
retry. Current three-cell figures do not imply that live or missing jobs scored.

The figure contains no timing comparisons. Any resource tables using these
experiments must identify observed shared-Threadripper timing under competing
CPU, memory-bandwidth and I/O demand, with unknown, potentially tool-dependent
effects. No isolated timing, recovered-timing repair, selection of fastest
repeats, general OrthoFinder superiority or publication readiness is claimed.
