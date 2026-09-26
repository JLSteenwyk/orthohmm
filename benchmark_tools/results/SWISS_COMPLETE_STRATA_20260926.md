# Complete Corrected SwissTrees Descriptive Strata

The admitted eight-method count report
[22178](qfo_recovered_swiss_uncertainty_22178.json), SHA256
`a599749d66433ec211ed1ab0e3a6a2a75eb99cb0dc57abc47c2d267930cbbe34`,
now supplies the previously missing recovered OrthoMCL rows in four exports:

- [Composition and sequence-property bins](swiss_descriptive_strata_20260926/scores.md): 88 rows.
- [Alignment-identity bins](swiss_identity_strata_20260926/scores.md): 32 rows.
- [Historical fragment-annotation bins](swiss_fragment_strata_20260926/scores.md): 56 rows.
- [Mapped-tree duplication bins](swiss_duplication_strata_20260926/scores.md): 32 rows.

Each directory includes a machine-readable manifest and TSV. All bins,
feature admissions and statistic definitions remain unchanged. All seven
previously admitted methods' complete row dictionaries match the September
23 exports exactly. Empty bins remain missing, not zero.

OrthoMCL F1 is 0.689327/0.815408 in lower/higher entropy bins,
0.672216/0.834263 in lower/higher identity bins, and 0.802575/0.728505
in lower/upper mapped duplication-fraction bins. Historical annotation-positive
families have F1 0.758627 versus 0.769521 for fully matched unflagged families.
These are descriptive macro-precision/recall harmonic means, not averages
of family F1, independent validation, causal explanations or significance
tests. Unflagged proteins are not proven complete, identity is not calibrated
evolutionary distance, and mapped-tree annotations are reference-dependent.

## Reproduction

The four existing `benchmark_tools.export_swiss_*_strata` modules now accept
`--counts-sha256`; omitting it preserves the historical seven-method digest.
Use the same feature/admission arguments recorded in each manifest, the
22178 counts above, its explicit digest, and a fresh output directory.
For example, from the repository root:

```bash
python -m benchmark_tools.export_swiss_descriptive_strata \
  --counts benchmark_tools/results/qfo_recovered_swiss_uncertainty_22178.json \
  --counts-sha256 a599749d66433ec211ed1ab0e3a6a2a75eb99cb0dc57abc47c2d267930cbbe34 \
  --strata benchmark_tools/results/corrected_swiss_sequence_strata_20260918.json \
  --output /tmp/swiss-descriptive-complete
```

The recorded feature/source checks remain enforced. The helper digest changes
only to account for explicit count selection; scoring and binning code is
unchanged. Twenty-nine focused tests passed, including all four new exports,
exact prior-row equality, unchanged default digests and wrong-hash rejection.
This reproduces descriptive tables from admitted counts, not native inference.
Figures and manuscript strata prose still need integration.
