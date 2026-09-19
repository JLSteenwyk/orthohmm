# Corrected QfO Figure Supplement

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
