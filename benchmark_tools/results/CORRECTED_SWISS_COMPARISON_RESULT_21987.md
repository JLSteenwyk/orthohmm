# Corrected SwissTrees Comparator Intervals

## Provenance and Execution

Job 21987 completed exit 0:0 in 23 seconds with 2 CPUs and 64 GiB requested
on bizon. The two-hour/no-requeue job used frozen, pushed executor
`10338e2a5046e522f2c1990e791532d3376e8982`, retained at
`benchmarks/work/publication_corrected_swiss_comparison_v1`. Terminal
accounting preceded result inspection; the batch log is empty. Runner SHA-256:
`d757c246c1b6311a9fb5a885780d2fc9f0dbbc167448c8e49019a4ff7915aad9`.

The [execution note](CORRECTED_SWISS_COMPARISON_EXECUTION_20260919.md)
preserves both original protocol hashes and the six-admitted-method
comparison-manifest hash. All eight methods and contrasts remain represented.
Corrected raw family counts were reconstructed against the same 18-family,
10,765-relation reference universe; old prediction counts were not reused.
The analysis uses 100,000 shared PCG64 family draws, seed 20260920, and
Bonferroni adjustment over all 24 planned endpoints, including unavailable
comparisons. Native macro precision/recall and harmonic F1 are recomputed
within each draw. No tuning or new contrast was introduced.

Submitted command:

```bash
sbatch --parsable --job-name=qfo_corrected_comparator_uncertainty \
  --nodelist=bizon --cpus-per-task=2 --mem=64G --time=02:00:00 --no-requeue \
  --output=benchmarks/work/qfo_corrected_comparator_uncertainty_%j.log \
  --wrap='env -u PYTHONPATH -u PYTHONHOME -u LD_PRELOAD -u LD_LIBRARY_PATH PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/bizon/anaconda3/bin/python -B benchmarks/work/publication_corrected_swiss_comparison_v1/benchmark_tools/run_corrected_swiss_comparison.py --comparison benchmark_tools/results/qfo_corrected_comparison_20260919_v5/manifest.json --comparison-sha256 7256920214d2eeda9cf596c7c39212992c7f3fcc36dd61e1848fa44d5725e01e --baseline benchmark_tools/results/qfo_swiss_counts_20260917.json --protocol benchmark_tools/results/QFO_SWISS_COMPARATOR_UNCERTAINTY_PROTOCOL_20260917.md --release-protocol benchmark_tools/results/QFO_CORRECTED_RELEASE_PROTOCOL_20260917.md --output benchmarks/work/qfo_corrected_comparator_uncertainty_v1.json'
```

[Retained result](qfo_corrected_comparator_uncertainty_21987.json):
213,102 bytes, SHA-256
`b05b52a50a1361537f9d4aef581f8fb8366db1591cb3cadd7aff5bc074953c95`.
The scheduled output and retained copy match. Six contrasts/18 endpoints
are estimable; the two missing-method contrasts remain null. The report
does not declare the complete panel or publication package complete.

## Results

Candidate minus reference; values are raw 0-to-1 F1 differences, not percent
relative improvement. The intervals below use the full 24-endpoint adjustment.

| Candidate | Reference | F1 difference | Adjusted interval | Family wins/ties/losses |
| --- | --- | ---: | --- | --- |
| OrthoHMM high-sensitivity | OrthoFinder full | -0.162915 | [-0.270844, -0.020071] | 1/0/17 |
| OrthoHMM phylogenetic | OrthoFinder full | -0.014900 | [-0.087078, 0.072090] | 3/6/9 |
| OrthoFinder sequence-only | OrthoFinder full | -0.157898 | [-0.317472, 0.023692] | 2/1/15 |
| SonicParanoid | OrthoFinder full | -0.049954 | [-0.135157, 0.038111] | 3/1/14 |
| Proteinortho | OrthoFinder full | -0.130302 | [-0.255797, -0.020741] | 3/2/13 |
| FastOMA | OrthoFinder full | Not admitted | Not estimable | Not estimable |
| OrthoMCL | OrthoFinder full | Not admitted | Not estimable | Not estimable |
| OrthoHMM phylogenetic | OrthoHMM high-sensitivity | 0.148015 | [0.049791, 0.266015] | 16/1/1 |

All three adjusted intervals for phylogenetic OrthoHMM versus full OrthoFinder
include zero. This does not show superiority or equivalence. Relative to
high sensitivity, phylogenetic OrthoHMM improves F1 and precision under this
adjustment; its recall interval includes zero. Candidate expansion also
changes between these configurations, so this is not a pure reconciliation
effect. Full OrthoFinder has higher F1 than high-sensitivity OrthoHMM and
Proteinortho under this conditional analysis.

The OrthoFinder sequence-only checkpoint shows lower precision and higher
recall than its full pipeline, with both adjusted intervals excluding zero;
its adjusted F1 interval includes zero. SonicParanoid's three adjusted
intervals all include zero. Across the 18 available endpoints, eight adjusted
intervals exclude zero. All nominal intervals, adjusted intervals and family
differences are retained, including unfavorable and null results.

These are development-exposed, conditional intervals over only 18 curated
families. Shared evolutionary history and merged predictions can violate
exchangeability despite disjoint represented proteins. The adjustment does
not cover past method selection or all publication analyses. No intervals
are transferred to TreeFam, other QfO endpoints, the secondary mean or
historical inputs; a significance change between releases is not itself a
release-by-method interaction test.

## Reproduction and Figure

[Independent numerical reproduction](qfo_corrected_comparator_reproduction_21987.json)
recalculates smoothed raw-count probabilities and weighted family sums
without importing the production bootstrap/statistic helpers. All 18
available endpoints, points, both interval types, per-family differences
and wins/ties/losses agree within `1e-12`; missing rows and panel completeness
are also checked. The same NumPy RNG/quantile implementation is used. This
does not rerun inference, scoring or independently establish raw truth.
Reproduction/plotting source is committed and pushed at `87e7a49`.

The [PNG](corrected_swiss_comparison_figure_21987/corrected_swiss_comparison.png),
[PDF](corrected_swiss_comparison_figure_21987/corrected_swiss_comparison.pdf),
and [SVG](corrected_swiss_comparison_figure_21987/corrected_swiss_comparison.svg)
show all eight contrasts across three metrics. Thin intervals are adjusted;
thick intervals are nominal. Plot values alone are multiplied by 100 for
percentage-point display. The [TSV](corrected_swiss_comparison_figure_21987/endpoints.tsv)
and [full table](corrected_swiss_comparison_figure_21987/endpoints.md) retain
raw units and all 24 planned endpoint rows, including six null rows.
The figure manifest binds source/input/output hashes. The rendered PNG was
visually inspected: intervals fit the axes and labels are not clipped or
overlapping. Pixel checks also confirm a nonblank 2700-by-1260 image.

Tests: 84 audit/kernel/export cases and 48 kernel/driver/reproduction/plot
cases passed (overlapping suites, not additive). Synthetic full-panel
outputs exactly match the frozen historical algorithm when given identical
counts. Failure tests cover changed input/protocol/source, count identities,
missing comparators, multiplicity, incorrect arithmetic and output overwrite.
Generated Matplotlib SVG retains its original serialization, including
standard trailing whitespace, so its recorded hash remains intact.
The TSV also retains trailing tab delimiters for its six missing-result
rows: those delimiters encode empty columns, not zero scores. Git whitespace
checks flag those generated rows; stripping them would change the table
schema. Manual Markdown, Python and JSON whitespace checks pass separately.

FastOMA and OrthoMCL admissions, the final eight-method refresh, secondary
strata, and consolidation into the complete publication bundle remain open.
