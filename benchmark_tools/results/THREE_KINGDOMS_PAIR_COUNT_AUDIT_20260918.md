# Independent Three Kingdoms Pair-Count Audit

All eight retained normalized-group results reproduce their historical
summary counts. The new auditor does not import the scorer, its group reader,
or the summary parser. It builds reference/prediction partition intersections:
TP is the sum of choose-two intersection sizes, predicted positives are the
sum of choose-two reference-restricted predicted-group sizes, and truth is
the sum of choose-two reference-group sizes. FP and FN follow by subtraction.
This is independent of the scorer's explicit pair enumeration.

| Method | TP | FP | FN | F1 |
| --- | ---: | ---: | ---: | ---: |
| OrthoHMM high sensitivity | 5176 | 0 | 2176 | 0.826309 |
| OrthoHMM satellite phylogeny | 5685 | 0 | 1667 | 0.872133 |
| OrthoFinder 3.1.5 full | 7186 | 0 | 166 | 0.988582 |
| OrthoFinder sequence-only diagnostic | 7269 | 72 | 83 | 0.989451 |
| ProteinOrtho 6.3.6 | 6396 | 0 | 956 | 0.930463 |
| SonicParanoid 2.0.9 historical input | 7265 | 48 | 87 | 0.990794 |
| FastOMA 0.3.5 | 5565 | 0 | 1787 | 0.861655 |
| OrthoMCL 1.4 | 7191 | 0 | 161 | 0.988929 |

Every method has 7352 reference pairs. Counts and coverage match exactly;
full-precision metrics match within 1e-12. Historical group, score-file and
reference hashes match the panel; all evidence is rehashed after calculation.
The machine-readable result is `three_kingdoms_pair_counts_audit_20260918.json`.

Nineteen focused tests pass, including 100 seeded random partitions checked
against explicit pair sets, missing genes, merges, splits, singleton groups,
out-of-reference genes, malformed partitions, nonfinite scores, source
mutations and corrupted reported counts.

This verifies normalized-group arithmetic against the historical summary,
not native conversion, historical input consumption, or runtime provenance.
The historical SonicParanoid input mismatch remains unresolved until the
matched rerun is validated. Out-of-reference predictions are not penalized;
the metric includes within-species co-membership pairs and is not a
proteome-wide pairwise-orthology test. No new method accuracy admission or
table replacement follows from this audit alone.

Reproduce from repository root:

```bash
python benchmark_tools/audit_three_kingdoms_pair_counts.py --repo . \
  --panel benchmark_tools/results/three_kingdoms_parity_20260907.json \
  --reference three_kingdoms/busco/reference_orthogroups.txt \
  --output /tmp/three_kingdoms_pair_counts_audit.json
```

The output destination must not already exist.
