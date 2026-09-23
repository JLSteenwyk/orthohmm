# Corrected QfO Figure Supplement

## Seven-Method Comparator Supplement (2026-09-23)

The separate FastOMA-inclusive comparison supplement is frozen at
`27113f40df3c9f21d4e9e294ac1bc079c2b2ba8d`. Earlier bundles below remain
unchanged. It includes seven admitted methods and the explicitly missing
OrthoMCL contrast: 21 estimated endpoints out of 24 planned. This supplement
does not replace or revalidate the factorial and composition-strata bundles.

The numerical helper has gained portable reproduction since the original
22098 figure was rendered. A new export binds its current source identity;
the PNG, endpoint TSV and endpoint Markdown remain byte-identical to 22098.
The PNG was visually inspected. PDF/SVG metadata differ after regeneration.
Fifty-four focused audit, bundle and plotter tests passed. Whitespace checks
pass for source, JSON and Markdown; generated SVG path whitespace and empty
TSV cells for unavailable OrthoMCL values are preserved without rewriting.

- [Direct-byte audit](publication_fastoma_figure_integrity_20260923.json):
  all recorded identities match, no untracked dependencies.
- [Bundle manifest](publication_fastoma_figure_bundle_20260923.json):
  SHA-256 `8f284a1e5387cd107efbfc53bf9614a303657cc417106de236f39d4c8a8d8120`.
- [Relocated verification](publication_fastoma_figure_bundle_relocation_20260923.json):
  12 files, five outputs, 609,328 bytes excluding the 7,282-byte manifest.
- [Bundled arithmetic reproduction](publication_fastoma_bundle_arithmetic_20260923.json):
  all 21 estimated endpoints reproduced within 1e-12 using Python 3.10.13
  and NumPy 2.2.6. Historical absolute paths were not accessed.
- Local archive: benchmarks/work/publication_fastoma_figure_evidence_20260923_27113f4.tar.gz,
  220,341 bytes, SHA-256
  `f8f23387a82281350f344ee3b79304bbc7d98bd07d2e6a08f1ddce919047d467`.

Extracted to `/tmp/orthohmm-fastoma-figure-PmU2F5Ho` and checked there with
`/usr/bin/python3 -I benchmark_tools/bundle_publication_figures.py verify .`.
The bundled numerical verifier also ran from that directory with `-I -B`,
using the already pinned Swiss analysis environment. Its output was written
outside the bundle so the exact file inventory remains intact.

Reproduce with a fresh output path:

```sh
python -m benchmark_tools.bundle_publication_figures build \
  --repo . --revision 27113f40df3c9f21d4e9e294ac1bc079c2b2ba8d \
  --audit benchmark_tools/results/publication_fastoma_figure_integrity_20260923.json \
  --output /tmp/fastoma-figure-evidence
python -I /tmp/fastoma-figure-evidence/benchmark_tools/bundle_publication_figures.py \
  verify /tmp/fastoma-figure-evidence
python -I -B /tmp/fastoma-figure-evidence/benchmark_tools/reproduce_corrected_swiss_comparison.py \
  --retained-counts-only \
  --results /tmp/fastoma-figure-evidence/benchmark_tools/results/qfo_fastoma_swiss_uncertainty_22098.json \
  --results-sha256 121cc8adf3cd63878f19006ec5500a13f879d042eccd565e1ae5f5a863f43fb1 \
  --output /tmp/fastoma-figure-arithmetic.json
```

The final command requires NumPy; use the pinned environment documented in
the reproduction guide. This is direct evidence preservation and retained
count arithmetic, not raw scoring, inference, plot regeneration outside the
checkout, cross-platform validation, full dependency packaging, cleared
redistribution, external deposition or publication readiness. DGX timing is
deferred by user request; no timing claim follows from this supplement.

## Earlier Three-Figure Export

The expanded supplement, frozen at commit
`5132ddaee322ce51c19c7ac7cb0fe7cb172e0e56`, contains the corrected factorial,
composition-stratified SwissTrees, and six-admitted-method SwissTrees
comparison figures. The original single-figure export below is preserved.
There are 27 files, 12 outputs and 1,919,243 bytes excluding the 16,022-byte
bundle manifest. Missing comparator results remain explicitly unavailable.

- [Direct-byte audit](publication_corrected_figure_integrity_20260919_v2.json):
  SHA-256 `3ba8464c42526c513f8e578551d4a64d3ec3a318745f6cc6599558ae7aebc7dd`.
- [Bundle manifest](publication_corrected_figure_bundle_20260919_v2.json):
  SHA-256 `a71a3bf4dbdcb13d70110d6dd076ef303791b81bc0db3de700d57c6b2ea25cd2`.
- [Relocated verification](publication_corrected_figure_bundle_relocation_20260919_v2.json):
  SHA-256 `42c48bd4cd67b6cadca401cb83340eb4e76ec117e40947a626b2f4058e6a8a50`.
- Local archive:
  `benchmarks/work/publication_corrected_figure_evidence_20260919_5132dda.tar.gz`,
  804,449 bytes, SHA-256
  `4f64f9bebfe3618a322f501c27f873356a65217c06c79ced6955c76a9c5a8f6b`.

The archive was extracted to `/tmp/orthohmm-corrected-figures-v2.H6fguz`.
The bundled verifier passed with `/usr/bin/python3 -I`, using that directory
as its working directory and `.` as the verification target. It resolves
recorded dependencies within the export, not historical absolute paths.
The corrected-current audit scope checks all three figures; the historical
and original corrected-factorial scopes remain separate. Forty-five focused
audit and bundle tests passed when this scope was added.

Reproduce the expanded export using a fresh output path:

```sh
python -m benchmark_tools.bundle_publication_figures build \
  --repo . --revision 5132ddaee322ce51c19c7ac7cb0fe7cb172e0e56 \
  --audit benchmark_tools/results/publication_corrected_figure_integrity_20260919_v2.json \
  --output /tmp/corrected-three-figure-evidence
python -I /tmp/corrected-three-figure-evidence/benchmark_tools/bundle_publication_figures.py \
  verify /tmp/corrected-three-figure-evidence
```

This is a direct-file integrity and relocation check only. It does not
regenerate statistics or plots, rerun native tools, recursively package
all dependencies, or resolve redistribution rights. The archive remains
local; no external deposit or DOI has been created. Publication readiness
remains false.

## Original Single-Figure Export

The corrected-release factorial figure now has a separate relocatable
direct-evidence bundle, built from commit `706fc03`. The historical 18-panel
bundle remains unchanged. Separate bundles preserve the distinct historical
and corrected plotter versions without rewriting their recorded identities.

The supplement contains one four-panel figure in PNG/PDF/SVG, its original
manifest, recorded source result, plotter and two statistical helpers, plus
the committed audit, repository license and standalone verifier: 11 files,
626,751 bytes, excluding the 6,430-byte bundle manifest.

- [Direct-byte audit](publication_corrected_figure_integrity_20260919.json):
  SHA-256 `f46dab5347a63a3e21a43de0b707a6da9403a32a97990435e21e8a7363f6d4fa`.
- [Bundle manifest](publication_corrected_figure_bundle_20260919.json):
  SHA-256 `1b2678006272206324482cb6e07d7ebe6b424be7475d33ef432257a09c178e3f`.
- [Initial verification](publication_corrected_figure_bundle_verification_20260919.json)
  and [relocated verification](publication_corrected_figure_bundle_relocation_20260919.json)
  are byte-identical.
- Local archive:
  `benchmarks/work/publication_corrected_figure_evidence_20260919_706fc03.tar.gz`,
  321,319 bytes, SHA-256
  `5f861546e638f8996f046e001de8ea1f4b9d582fda84ac1854f87cc28fadfc9b`.

The archive was extracted outside the repository at
`/tmp/orthohmm-corrected-figure-XeB5fp`. Its bundled verifier passed with
`python -I` and working directory `/tmp`; verification reads explicit relative
relocations, not the source files at historical absolute paths. Fifty-five
focused audit, bundle and factorial-plot tests passed.

Reproduce using fresh output paths:

```sh
python -m benchmark_tools.bundle_publication_figures build \
  --repo . --revision 706fc03 \
  --audit benchmark_tools/results/publication_corrected_figure_integrity_20260919.json \
  --output /tmp/corrected-factorial-evidence
python -I /tmp/corrected-factorial-evidence/benchmark_tools/bundle_publication_figures.py \
  verify /tmp/corrected-factorial-evidence
```

This proves direct-byte preservation and relocation, not regenerated plots,
statistics, scoring or native inference. Transitive dependencies, raw data,
the evolving corrected comparator table and redistribution review remain
outside this supplement. No external deposition or DOI was created, and no
new accuracy, timing or publication-readiness claim follows.
