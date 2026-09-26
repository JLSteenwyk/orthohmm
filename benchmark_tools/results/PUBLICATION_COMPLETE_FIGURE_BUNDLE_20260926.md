# Complete SwissTrees Figure Evidence Bundle

The explicit `corrected-complete` figure audit selects the September 26
eight-method comparison, leaving historical inventories unchanged. It checks
the admitted source, complete-panel flag, eight estimates/contrasts and 24
multiplicity endpoints in addition to recorded file identities. Sixty-two
figure-audit and bundle tests pass, including partial-panel rejection cases.

The bundle was built from commit
`81645d2c2c0c523b6925fdc1529f72525bbde5fa` and copied to
`/tmp/orthohmm-complete-figure-ASEAtI/bundle`. Its standalone standard-library
verifier ran there with isolated Python (`-I -B`), outside the checkout.
[Build receipt](publication_complete_figure_bundle_20260926.json) and
[relocation receipt](publication_complete_figure_bundle_relocation_20260926.json)
agree: 12 files, 954,077 payload bytes, one panel and five outputs
(PNG/PDF/SVG plus Markdown/TSV endpoints). The bundle manifest SHA256 is
`f9c5a6d5cdcc71eb7bad9066b6eb078fd0f482a33dcd79bb3290b571e015b796`.

The bundled standalone statistical verifier also reproduced all 24 endpoints
within 1e-12 from the bundled count report using the existing clean Python
3.10.13/NumPy 2.2.6 environment. See the
[arithmetic receipt](publication_complete_bundle_arithmetic_20260926.json).

## Reproduce

From the repository root, choose a fresh output path:

```bash
python -m benchmark_tools.bundle_publication_figures build \
  --repo . --revision 81645d2c2c0c523b6925fdc1529f72525bbde5fa \
  --audit benchmark_tools/results/publication_complete_figure_integrity_20260926.json \
  --output /tmp/complete-swiss-bundle-new
python -I -B /tmp/complete-swiss-bundle-new/benchmark_tools/bundle_publication_figures.py \
  verify /tmp/complete-swiss-bundle-new
python -I -B /tmp/complete-swiss-bundle-new/benchmark_tools/reproduce_corrected_swiss_comparison.py \
  --results /tmp/complete-swiss-bundle-new/benchmark_tools/results/qfo_recovered_swiss_uncertainty_22178.json \
  --results-sha256 a599749d66433ec211ed1ab0e3a6a2a75eb99cb0dc57abc47c2d267930cbbe34 \
  --retained-counts-only --output /tmp/complete-swiss-arithmetic-new.json
```

NumPy is required for the final command, not for byte verification. Keep
new reports outside the bundle: verification rejects extra files.

This is direct figure evidence and retained-count arithmetic, not native
inference, raw-score admission, plotting regeneration, cross-platform
validation, third-party rights clearance or the full publication archive.
The full goal remains incomplete, including controlled resource comparisons.
