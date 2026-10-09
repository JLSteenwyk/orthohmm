# Terminal-Failure Manuscript Review

The new reporting revision retains four admitted native QfO cells and three
missing cells. It does not change the earlier scores, uncertainty or figures.
Native11 completed inference/conversion but scoring failed OUT_OF_MEMORY;
native12 inference exited -11 (SIGSEGV). Partial native12 phylogeny outputs
were not converted or scored. No inference retry was performed.

## Evidence

- [Reporting snapshot](native_qfo_terminal_failures_20261009_v1/report.json),
  SHA256 `863ebc334dbdb5a41d1a9e9bba56be82896f0ff26c98366047928c7ef802966f`.
- [Manuscript source](native_qfo_terminal_failures_20261009_v1_manuscript.md),
  [HTML review](native_qfo_terminal_failures_20261009_v1_review.html), and
  [PDF review](native_qfo_terminal_failures_20261009_v1_print/document.pdf).
  PDF SHA256 `470638f386ef64af6ee44061c22410dabeda40db93b5ff36c50068bb8575c177`.
- [Render/asset receipt](native_qfo_terminal_failures_20261009_v1_html_assets.json),
  [print receipt](native_qfo_terminal_failures_20261009_v1_print/print.json),
  [page-bounds receipt](native_qfo_terminal_failures_20261009_v1_layout_review/report.json),
  and [citation inventory](native_qfo_terminal_failures_20261009_v1_citations.json).
- [Independent content readback](native_qfo_terminal_failures_20261009_v1_content_readback.json),
  SHA256 `c27e5597b173775db639e4a42359359d8ff88e9de6d4858c3242927fd7eaaa25`.
  It retains the executed standard-library/PyMuPDF validation code.

## Actual Checks

The independent JSON comparison found all four admitted rows identical to
the frozen snapshot and all missing endpoint scores/means null. CSV readback
checked 42 endpoint rows and seven status rows, with 18 missing values.
Removing only the new terminal-failure section reproduces the original
parent manuscript body exactly.

The actual new PDF has 22 pages. All 22 page PNGs were visually inspected in
this continuation, including the status table and failure descriptions on
pages 10-11. No clipped or overlapping text was observed. The automated bounds
check found zero violations; that check alone would not prove visual fidelity.
The evidence/hash paragraph is compact, and final journal-specific typesetting
remains outstanding. This is not a submission-ready layout certification.

Independent PDF text checks found all 24 admitted endpoint values and eight
selected failure/resource phrases. All ten directly linked PDF assets were
freshly decoded on every page and checked nonblank. This is asset decoding,
not a new manual visual assessment of every historical figure or replotting.

The render checked 114 local link occurrences and 110 unique targets. All
19 used citation ids resolved against the October bibliography, including the
corrected IQ-TREE 3 reference; the separate inventory found no unresolved
explicit citations. Resolution does not prove attribution adequacy or rights.
The 36,151,387-byte native12 review manifest remains an untracked direct target
at this checkpoint; no portable direct-review archive for this revision has
yet been built or restored.

## Limits And Next Work

The source header records its pre-render workflow state. These subsequent
dated receipts establish this local render/readback only; they do not alter
the frozen source or imply whole-study release/deposition.

Root cause and crash location of native12 remain unknown. Missing native cells
prevent a complete fresh factorial. Development exposure, family-overlap,
functional/TreeFam uncertainty, raw provenance, runtime restoration and
redistribution limitations remain explicit in the manuscript.

Shared-host timings remain potentially confounded observations, not isolated
efficiency rankings. No new score, bootstrap draw, default, independent
confirmation, raw scientific admission, redistribution clearance or archival
DOI follows from this review. Next integrate the new dated review into a
truthfully scoped reporting component and verify its actual copied payloads.
