# Native QfO Figure Result

The real export after the validated repair in commit `50fdea48` produces
[PNG](native_qfo_p0c0_figure_20261006_v1/native_qfo_p0c0.png),
[vector PDF](native_qfo_p0c0_figure_20261006_v1/native_qfo_p0c0.pdf) and
[SVG](native_qfo_p0c0_figure_20261006_v1/native_qfo_p0c0.svg), plus full-precision
[12-point score table](native_qfo_p0c0_figure_20261006_v1/scores.tsv) and
[three-effect interval table](native_qfo_p0c0_figure_20261006_v1/swiss_intervals.tsv).
The [manifest](native_qfo_p0c0_figure_20261006_v1/manifest.json) is 36,793 bytes,
SHA256 `6d21105f36311f2529ee9447d93790c40125f34216be6b0e8daea69aa3180577`.
It records exact original-environment binding replay, separate rendering
versions, all 114 evidence records and five output hashes. No inference,
scoring, count admission, bootstrap or native timing is repeated.

## Content And Visual Review

The independent [content readback](native_qfo_p0c0_figure_readback_20261006.json)
checks all score and interval TSV values exactly against the original snapshot
and binding, rather than using the plotting helper's data function. It also
checks source/evidence/output hashes, SVG labels, PDF header and image pixels.
PNG dimensions are 2500 by 1600 with 157,392 nonwhite pixels, 1,972 exact teal
pixels and 3,311 exact brick pixels. The PDF raster preview is 1600 by 1024
with 69,090 nonwhite pixels and both method colors.

Both the original PNG and an actual `pdftoppm` PDF rasterization were inspected
using `view_image` on 6 October 2026. All four panels, legends, axes and three
limitation lines are visible without clipping or incoherent overlap. The
small TreeFam-A difference makes its two markers overlap at this scale; the
exact tables and manuscript explicitly retain the decrease. Panel D visibly
includes zero in the adjusted F1 interval. This separate inspection does not
rewrite the original manifest's `visual_review_complete: false` or turn
nonblank pixels into automatic visual certification.

The [retained PDF preview](native_qfo_p0c0_pdf_preview_20261006.png) is byte-identical
to the temporary preview referenced in the content receipt (112,953 bytes,
SHA256 `85ad17e9999b4f69a91fccfa0a386bed2264dda48013755dca44e534b0069791`).
The generated preview is copied without overwriting an existing destination;
the copy command emits a non-portability warning about `cp -n`, not an analysis
failure. Recreate previews with a fresh destination using:

```bash
pdftoppm -singlefile -scale-to 1600 -png \
  benchmark_tools/results/native_qfo_p0c0_figure_20261006_v1/native_qfo_p0c0.pdf \
  /tmp/orthohmm_native_qfo_p0c0_pdf_new_review
benchmarks/work/release_alert_refresh_20261001/venv/bin/python -B \
  -m benchmark_tools.review_native_qfo_figure \
  --figure-directory benchmark_tools/results/native_qfo_p0c0_figure_20261006_v1 \
  --pdf-preview /tmp/orthohmm_native_qfo_p0c0_pdf_new_review.png \
  --output /tmp/orthohmm_native_qfo_new_content_review.json
```

## Interpretation And Limits

Two of seven fresh native QfO cells are admitted; five remain unavailable.
Initial HMM search is on in both P0/C0 cells. R changes group-clique predictions
to inferred pairs, so this is neither a non-HMM control nor a selected-default
OrthoHMM-versus-OrthoFinder comparison. Reconciliation increases the observed
SwissTrees F1 by 10.0390 percentage points, but the adjusted interval
[-4.6418, 24.4199] includes zero. Precision rises while recall falls. Only one
contrast has exact matched family records allowing reuse of the retained
100,000 draws; all 42 endpoints remain in adjustment. These 18 families are
development-exposed, with conditional exchangeability/percentile limits.

TreeFam-A F1 decreases. No paired intervals for other QfO endpoints or the
secondary mean are added; unseeded FAS sampling and dependence remain. Original
22437 timing stays failed/null/ineligible, and this figure has no timing panel.
All reported timings elsewhere remain shared-Threadripper observations with
unknown, potentially tool-dependent CPU, memory-bandwidth and I/O contention.
The working manuscript/checklist/guide incorporate this bounded evidence;
the earlier rc4 archive/PDF payloads remain historical and unchanged.
The full publication goal remains active and readiness unproven.

## Verification

The joined suite passes 169 tests in 8.56 seconds. Seventeen new content-review
cases check exact artifacts, current manuscript/caption numbers, incomplete
scope and rejection of changed scores, statistic labels, duplicate/missing
rows, intervals, SVG labels, blank pixels, PDF header, hashes and admission
flags. Refreshed output hashes do not bypass the semantic checks. The actual
retained-environment binding replay also passes. These tests validate the
reporting workflow and retained artifacts, not scientific generalization or
the unfinished full-goal requirements.
