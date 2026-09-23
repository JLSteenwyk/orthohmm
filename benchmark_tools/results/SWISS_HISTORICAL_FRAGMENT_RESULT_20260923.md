# Historical Fragment Annotation Results

## Admission and Provenance

The [prespecified protocol](SWISS_HISTORICAL_FRAGMENT_PROTOCOL_20260923.md)
selected historical UniSave records without using prediction outcomes.
Collection job 22116 completed 0:0 in 24:16; independent admission job 22117
completed 0:0 in 00:12. All 563 proteins matched accession, taxon, sequence
version and exact sequence digest. Of these, 549 matched the baseline-release
interval and 14 required later records for the input sequence version.
There are 277 Reviewed and 286 Unreviewed entries. Annotation admission is
not experimental validation of biological completeness.

The [admission report](swiss_historical_fragment_admission_22117.json) has
SHA256 `a480f32666bb96de3395213c7cad024c170318ca110c389f15e4756f8052adf8`.
Raw histories, entries and acquisition receipts remain under
`benchmarks/work/swiss_historical_fragment_panel_20260923` and are fingerprinted
in that report. The collector and verifier use separate selection logic but
share Bio.SwissProt parsing; this is not a wholly independent parser audit.

Eleven proteins carry positive fragment/incomplete-sequence annotations:
A7RPH5, B7PU93, B8BVL0, C3XR12, C3XTC7, C3XTC8, C3XTD0, C3XTD1, C3XTD2,
C3YHW0 and Q98924. The five annotation-positive families are APP, BAR, HOX,
NOX and PSEN. Thirteen families have all members matched and unflagged;
none has missing annotations without a positive member. Unflagged does not
mean proven complete. A family-level positive does not label every member
or every pair as fragmented.

## Descriptive Results

[Full generated tables](swiss_fragment_strata_20260923/scores.md),
[machine-readable values](swiss_fragment_strata_20260923/scores.tsv) and
[export manifest](swiss_fragment_strata_20260923/manifest.json) retain all
eight methods across seven views (56 rows), including missing OrthoMCL
and empty bins. F1 is the harmonic mean of macro precision and macro recall,
not the mean of family F1 scores. Values below are percentages.

| Method | Positive (5) | Unflagged (13) | Baseline-only unflagged (11) |
|---|---:|---:|---:|
| OrthoHMM high-sensitivity | 63.597 | 70.434 | 68.325 |
| OrthoHMM phylogenetic | 83.827 | 83.143 | 84.562 |
| OrthoFinder 3.1.5 full | 87.926 | 83.498 | 83.539 |
| OrthoFinder sequence-only | 64.410 | 70.295 | 66.133 |
| SonicParanoid | 81.459 | 79.148 | 78.457 |
| ProteinOrtho | 64.785 | 74.193 | 74.308 |
| FastOMA (supplied tree) | 79.794 | 77.323 | 76.827 |
| OrthoMCL 1.4 | NA | NA | NA |

Phylogenetic OrthoHMM minus full OrthoFinder is -4.099 percentage points
in positive families and -0.355 in unflagged families. In the positive bin,
OrthoHMM precision/recall is 98.736%/72.829%, versus 91.377%/84.726% for
OrthoFinder: the higher precision does not offset lower recall.

The baseline-only sensitivity analysis treats all 14 later-version entries
as missing. Positive membership stays unchanged; ASTER and VATB move from
unflagged to missing-without-positive. The remaining eleven unflagged
families give a +1.023-point OrthoHMM difference, while the two now-missing
families give -9.476 points. This sign change reflects a different family
composition, not an improvement of the method or evidence of superiority.
The overall 18-family values remain 83.351% versus 84.841%.

## Scope and Reproduction

These development-exposed associations have no new confidence intervals,
significance tests, causal claims or independent generalization claims.
Family size, divergence, domains and taxonomic membership may confound
the contrasts. Sequence-only OrthoFinder uses group-derived cliques, not
native phylogenetic ortholog predictions. FastOMA uses a supplied tree.
The earlier empty FASTA-header fragment bin is a different annotation
source, not a contradiction or a negative biological finding.

After restoring the fingerprinted local inputs, use a new output directory:

```bash
python -B -m benchmark_tools.export_swiss_fragment_strata \
  --counts benchmark_tools/results/qfo_fastoma_swiss_uncertainty_22098.json \
  --admission benchmark_tools/results/swiss_historical_fragment_admission_22117.json \
  --admission-sha a480f32666bb96de3395213c7cad024c170318ca110c389f15e4756f8052adf8 \
  --output benchmarks/work/swiss_fragment_strata_reproduction
```

The exporter checks the admitted bins and pinned statistic implementation,
then checks source records before and after export. It does not rerun
inference or modify the frozen scientific configuration. The dated rendered
manuscript remains a historical snapshot and does not yet include this result.
