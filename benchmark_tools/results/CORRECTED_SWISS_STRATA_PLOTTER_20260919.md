# Corrected SwissTrees Strata Figure Preparation

Added `benchmark_tools/plot_corrected_swiss_strata.py` for the result schema
produced by queued21896. No real corrected-strata result has been plotted.
The historical domain-strata plotter contains different schemas, bin labels
and outcome-specific text and is left unchanged.

The new plot requires the corrected result status, scientific-input and
uncertainty admission flags, frozen replicate/seed/multiplicity metadata,
18nonoverlapping families, the three planned contrasts and three metrics.
It validates interval eligibility, point availability, finite mathematical
bounds and nominal/adjusted interval ordering. It does not reconstruct raw
counts or independently establish the truth of an admission flag; those
checks remain upstream and result review is required before use.

All27primary endpoints are displayed: lower and higher composition entropy
bins plus their interaction, for three method contrasts and F1/precision/
recall. Axes expand symmetrically to include every point and adjusted bound;
there is no outcome-based endpoint selection or clipping. Unavailable
intervals are explicitly labeled and remain blank rather than zero in the
numeric table. The missing-composition bin remains descriptive, outside the
primary interval family, with its family count disclosed.

## Outputs and Interface

After terminal input validation and independent result review:

```sh
python -m benchmark_tools.plot_corrected_swiss_strata \
  --results REVIEWED_RESULT_JSON --results-sha256 REVIEWED_SHA256 \
  --output FRESH_FIGURE_DIRECTORY
```

The output directory must be fresh. It contains PNG/PDF/SVG figures,
`endpoints.tsv` with all27raw-scale endpoint values/bounds, and a manifest
binding the input checksum, plot source, Matplotlib version and four output
checksums. Only the visual display multiplies values by100. The source is
checked again before the manifest is written; failures leave partial evidence
rather than overwriting a prior result. Neither rendering nor the manifest
claims publication readiness.

## Validation

61focused tests passed in2.85seconds across the new plotter, corrected
driver/bootstrap and historical domain plotter suites. The18new plot tests
cover inventory, metadata, admissions, overlaps, invalid bounds/intervals,
missing values, lossless table export, hashes, occupied output paths and
display-unit conversion. Test prediction data and admission metadata are
synthetic fixtures, not real100000-replicate scientific results.

A synthetic preview was rendered and visually inspected at
`benchmarks/work/corrected_swiss_strata_synthetic_layout_20260919.png`, with
an explicit SYNTHETIC TEST ONLY title. It is not committed or used as a
publication result. Labels, intervals, panels and footnotes fit without
overlap; final panel headings expand PPV/TPR to Precision/Recall.
Actual-result rendering, inspection, figure-bundle updates and manuscript
integration await21896and its independent review. No DGX access occurred.
