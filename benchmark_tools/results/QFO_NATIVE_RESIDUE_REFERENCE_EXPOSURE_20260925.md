# Native Residue Reference Exposure

The seven reviewed single-O deletions were mapped injectively to the retained
QfO 2020 reference using the benchmark mapping and native Darwin structures.
This is a descriptive exposure audit, not search admission or an accuracy
counterfactual. These proteins are not labeled as failed queries.

| Reference evidence | Proteins exposed / seven |
| --- | ---: |
| SwissTrees mapped membership, 18 eligible cases | 0 |
| TreeFam-A mapped membership, one pooled retained case | 0 |
| VGNC incident reference pairs, 23,934 total | 0 |
| EC annotation | 7 |
| Experimental GO annotation under the benchmark evidence-code filter | 0 |
| FAS annotation entry with at least one feature | 7 |

All seven belong to `UP000002487_188937`. Each has one EC annotation
and one Pfam feature type. Accessions are P58865, P58866, Q8TN68, Q8TS72,
Q8TS73, Q8TTA5 and Q8TTA9. Protein numbers and feature counts are retained
in the machine-readable report.

Absence from a reference does not establish absence of an indirect grouping
effect. Conversely, an EC annotation or FAS feature does not establish that
the protein belongs to a scored predicted pair. Search-hit and final-group
tracing remain outstanding; no claim of negligible score impact is justified.
The pooled TreeFam reference does not recover original family-level units.

## Reproduction

From the repository root:

```bash
/home/bizon/anaconda3/bin/python -B -m benchmark_tools.audit_native_residue_reference_exposure \
  --root "$PWD" \
  --output "$PWD/benchmarks/work/qfo_native_residue_reference_exposure_20260925"
```

The output directory must not already exist. The completed command returned
zero; Darwin emitted its seven-target completion marker and exited zero.
All recorded reference, review and helper identities were checked after the
read-only analysis as well as before consumption. The source script and
native transcript identities are included in the report.

Retained report: `qfo_native_residue_reference_exposure_20260925.json`,
17,285 bytes, SHA256
`db1db3070eb9886d4a831ea58a057b2c032b4664c875feb67a9c7db72300bfd8`.
All 39 focused selection, native report parsing and reviewed-representation
tests passed. No original input, database or prediction was changed.
