# Missing SwissTrees Accessions In Retained Source FASTAs

The staged Xenopus tropicalis FASTA is byte-identical to the retained
extracted canonical FASTA: SHA-256
`924f7d0ab9f197d8e1108394ae82d8c947c08172890f6ddf91d17513711ad6e3`.
Both contain9,639sequences. The retained additional FASTA contains36,072.
Inspection of `qfo_benchmark/prepare_input.sh` shows that its staging rule
copies canonical FASTAs and explicitly skips additional/DNA files. The script
was inspected, not rerun; no input files were replaced.

## Located Identities

Eleven of the14SwissTrees reference accessions missing from frozen inputs
occur exactly in the additional FASTA. Each has a literal "Isoform of"
header naming a protein present in the canonical input. For nine, the named
protein has a different numeric identity in the frozen QfO mapping; two named
proteins are unmapped. None restores the missing reference ID by aliasing.

| Missing Reference Accession | Header-Named Canonical Protein | Canonical Mapping |
| --- | --- | --- |
| F6PXE7 | A0JM23 | Different numeric ID |
| F6QLB4 | A0A6I8RCX7 | Different numeric ID |
| F6SX31 | Q0IIS3 | Different numeric ID |
| F6YNC0 | A0A6I8Q122 | Different numeric ID |
| F7AMC8 | A0A6I8QW21 | Different numeric ID |
| F7BEI0 | F6QJZ3 | Different numeric ID |
| F7CLH6 | Q0IJ12 | Different numeric ID |
| F7D5A4 | A0A6I8SZP6 | Unmapped |
| F7D5I9 | B0BM40 | Different numeric ID |
| F7DF67 | A0A6I8Q6R4 | Different numeric ID |
| Q505H9 | F7BAU5 | Unmapped |

A0A6I8Q068, B1H1F6 and Q0VGW2 occur in neither retained FASTA. This audit
does not claim they are absent from every archive sidecar or external source.
The JSON retains both selected and named-target descriptions, lengths,
sequence hashes and numeric IDs without copying raw sequences into git.

## Interpretation And Reproduction

This localizes11missing exact accessions to additional-file exclusion in the
retained source layout; staging did not delete records from the retained
canonical FASTA. It does not yet authenticate those extracted bytes against
the original archive or establish what input representation QfO recommends
for this resource version. Investigate archive integrity and representative
selection before deciding on any separately frozen rerun or sensitivity study.

"Isoform of" is source annotation, not independent orthology truth. Some
descriptions/gene labels warrant caution; do not infer equivalence or silently
remap reference relations to named targets. The official scores and full
reference denominators remain unchanged. Historical comparator input parity
and broader input/reference coverage remain separate audit requirements.

```sh
python benchmark_tools/audit_swiss_additional_sequences.py \
  --aliases benchmark_tools/results/swiss_sequence_alias_audit_20260917.json \
  --canonical qfo_benchmark/proteomes/extracted/Eukaryota/UP000008143_8364.fasta \
  --additional qfo_benchmark/proteomes/extracted/Eukaryota/UP000008143_8364_additional.fasta \
  --staged qfo_benchmark/input/UP000008143_8364.fasta \
  --mapping qfo_benchmark/benchmark-webservice/reference_data/2020/mapping.json.gz \
  --output benchmark_tools/results/swiss_additional_sequence_audit_20260917.json
```

Report SHA-256:
`04cd54bf2b4be6d47afbe0b2f7d02ca2a492c0cbd1fe6c2ffed4fd1426cf8177`.
Tests cover selected parsing, exact header placement, missing accessions,
duplicate/empty/malformed entries, distinct numeric identities and mismatched
canonical versus staged files.
