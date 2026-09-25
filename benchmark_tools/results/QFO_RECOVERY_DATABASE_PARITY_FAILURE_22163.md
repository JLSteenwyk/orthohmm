# Recovered Search Database Parity Failure

Search admission 22163 failed 1:0 after 3:15 on 25 September 2026.
The exact-input gate correctly rejected the native database; no BPO or
inference was released. The merged candidate remains preserved.

## Evidence

Under `benchmarks/results/qfo_blast_replacement_search_admission_v1`:

- `report.json`: SHA256 `bae3f65eff58c5a0c72f4460e989f0879230e7a9932611d6c1eab8ec4f0b1886`.
- `database/report.json`: SHA256 `21ac4a1ed5504ab93de3511d190d5cf3b78cadb0ca8898264968d297b9e6faeb`.

The full extracted database retains all 984,137 records and their order.
984,130 sequences match exactly. Seven sequences each lose one `O` residue,
reducing residues from 440,246,934 to 440,246,927. Independently comparing the
seven source/dump sequence pairs confirmed that removing `O` gives the exact
dump sequence in every case; there are no other differences in those pairs.

| Sequence | Zero-Based Ordinal | Deleted Position (One-Based) |
| --- | ---: | ---: |
| sp\|P58865\|MTMB1_METAC | 583643 | 202 |
| sp\|P58866\|MTMB2_METAC | 583644 | 202 |
| sp\|Q8TN68\|MTBB3_METAC | 585798 | 356 |
| sp\|Q8TS72\|MTBB2_METAC | 587199 | 356 |
| sp\|Q8TS73\|MTTB2_METAC | 587200 | 334 |
| sp\|Q8TTA5\|MTBB1_METAC | 587581 | 356 |
| sp\|Q8TTA9\|MTTB1_METAC | 587585 | 334 |

The retained production `benchmarks/results/qfo_corrected_orthomcl_v1/work/formatdb.log`
reports seven illegal-character removals at the corresponding one-based
sequence ordinals. This is a native parser transformation, not exact parity.

## Isolated Native Probe

In fresh `/tmp/orthohmm_legacy_o_probe_20260925`, used the installed
`SOFTWARE/blast-2.2.13/bin/formatdb -i input.fa -p T -l formatdb.log`, then
`fastacmd -d input.fa -p T -D 1`. The fixture was:

```fasta
>with_o
ACDEFGHIKLMNOPQRSTVWY
>without_o
ACDEFGHIKLMNPQRSTVWY
>with_x
ACDEFGHIKLMNXPQRSTVWY
```

The formatter warned that one illegal `O` was removed. The first two dumped
sequences were identical (`ACDEFGHIKLMNPQRSTVWY`); the `X` control remained
unchanged. No production input, database or search output was modified.

## Required Follow-Up

Do not relabel this database as exactly matching the corrected FASTA. Before
allowing native-tool representation differences, implement and test an
explicit, narrowly checked transformation contract with the exact IDs,
positions, source/dump hashes, formatter evidence and unchanged-record count.
Retain the failed exact-parity admission and use a fresh validator output.
Unknown changes must continue to fail. Check the query parser separately;
database extraction alone does not establish query-side handling.

Assess reference exposure and available hits for the seven affected proteins
before claiming their impact is negligible. If native representation is
retained for comparator fidelity, label it explicitly and do not infer a
counterfactual corrected score. Alternatively, normalized-input sensitivity
runs must be separate and cannot silently replace the frozen comparator.
Whole-table HSP validation did not run after this gate failure and remains
required. No new production validation or inference job has been submitted.
