# Exact Classification Of Native Sequence Differences

Followed up all1,151differences from the original input/native database
sequence audit, rereading hash-checked source FASTAs and selected native
database entries. The audit preserves original hashes and makes no changes
to input sequences, identities or scores.

## Result

| Difference Class | Sequences | Location |
| --- | ---: | --- |
| Noncanonical symbols replaced only by X | 178 | 173outside XENTR;5in XENTR |
| Different lengths; not aligned | 973 | All XENTR |
| Other same-length differences | 0 | None |

Across the178same-length pairs, all315changed positions are B-to-X(46),
O-to-X(7), U-to-X(212) or Z-to-X(50). Every other position is identical.
Thus all observed non-XENTR sequence differences are fully accounted for
by these representation changes, not observed substitutions between standard
amino acids or length changes. The original byte differences remain real;
the data were not normalized into the exact-match count.

This demonstrates the exact input-to-database relationship. It does not
identify the historical preprocessing implementation or show that all
inference methods handle these symbols equivalently. No performance or
accuracy consequence is inferred from the315positions alone.

All973length-different sequences occur in Xenopus. They are deliberately
not zipped position-by-position or aligned by a new heuristic; their exact
source and native hashes/lengths remain recorded. Corrected-release validation
must resolve input identities and content separately. This supports treating
the Xenopus release issue distinctly from native residue representation.

## Reproduction And Validation

```sh
python benchmark_tools/classify_qfo_sequence_differences.py \
  --audit benchmark_tools/results/qfo_original_input_sequence_audit_20260917.json \
  --output benchmark_tools/results/qfo_sequence_difference_classes_20260917.json
```

Report SHA-256:
`451d7bd614e915a4ad338ad695ecda795d34b2d8676faf5556d0d739358ed2dc`.
Fourteen combined classification/sequence-reader tests pass. Tests ensure
standard-to-X, reverse-X, standard substitutions and other equal-length
changes are not mislabeled; unequal-length pairs are not truncated by zip.
The report retains every differing accession, numeric ID, species, source
hashes and per-sequence substitution counts without adding raw sequences.
