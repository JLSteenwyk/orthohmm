# Original Input Sequence Bytes Versus Native Scorer Database

Compared every mapped original input sequence with its numeric entry in the
retained `ServerIndexed.db`, pinned by the prior assessment environment to
SHA-256`f71e282f504d007306a34da98ab717c2c7fad0a28e4d503caa47a042a41ec9a8`.
The database contains984,137entries. Structured XML record parsing verifies
required fields; each mapped accession must occur in its expected native
entry's alias field. Source files are checked before and after reading.

## Results

| Category | Count |
| --- | ---: |
| Original input accessions | 976,504 |
| Mapped input sequences | 975,514 |
| Exact case-sensitive sequence matches | 974,363 |
| Sequence differences | 1,151 |
| Unmapped input accessions | 990 |
| Native entries without mapped input | 8,623 |

All8,623missing native entries are XENTR. Of1,151sequence differences,
978are XENTR and173are in23other species. There are973different-length
comparisons, all in XENTR, and178same-length comparisons, including all173
non-XENTR cases. The machine report retains every difference's accession,
numeric ID, species, source FASTA, lengths and sequence hashes.

## Interpretation

This extends the accession-coverage audit: even an accession recognized by
the scorer can carry a different original input sequence. Exact mismatch is
not itself proof of biological sequence replacement. In particular, check
native handling of nonstandard residues or sequence normalization before
interpreting the173same-length non-XENTR differences as release errors.
No such normalization was applied or assumed in this audit.

The intended corrected release still needs its own sequence-content check.
The numeric-coverage comparison queued as21688 does not establish byte
identity. Do not equate a fixed accession inventory with full input/scorer
compatibility. Original-release results remain preserved and explicitly
limited; no sequence, accession or benchmark score is changed here.

The parser reads native database records independently of Darwin. It checks
ordinal mapping against975,514input aliases and total984,137entries, but
does not validate GO/EC annotations, every external cross-reference, or all
native scoring algorithms. Sequence differences are reported, not discarded
or silently normalized into matches. Historical comparator input parity is
still a separate provenance question.

## Reproduction

```sh
python benchmark_tools/audit_qfo_input_sequences.py \
  --prepared benchmark_tools/results/qfo_factorial_prepared_20260917.json \
  --mapping qfo_benchmark/benchmark-webservice/reference_data/2020/mapping.json.gz \
  --database qfo_benchmark/benchmark-webservice/reference_data/2020/ServerIndexed.db \
  --database-sha256 f71e282f504d007306a34da98ab717c2c7fad0a28e4d503caa47a042a41ec9a8 \
  --output benchmark_tools/results/qfo_original_input_sequence_audit_20260917.json
```

Report SHA-256:
`480a09d0c60c95274d93fb49e15ae95b45a99d0cdcd17c72bf9b3af066cc12c4`.
Seven tests cover field structure, duplicate/missing fields and aliases,
whitespace rejection, case-sensitive sequences, exact/different/missing
comparisons, alias disagreement and incomplete reference counts.
