# Original QfO Input Versus Scorer Mapping Coverage

Verified report SHA-256:
`e13fe6d08c6908e9e17959929eb485de3ad5b8a95b40ade40cdfc66bd41da9dc`.

Applied the archive comparison workflow to the original archive, not the
corrected download. All78canonical archive FASTAs exactly match the frozen
input hashes. There are976,504unique input accessions and975,514unique
mapped numeric IDs. No duplicate canonical accessions or mapped IDs occurred.

The frozen scorer mapping contains984,137distinct positive numeric IDs.
Of these,8,623have no original canonical input accession. The native
one-based species-boundary convention places all8,623in XENTR, Xenopus
tropicalis. Conversely,990input accessions are unmapped; every one is in
`UP000008143_8364.fasta`. The other77proteomes have no unmapped input
accessions or missing reference numeric IDs under this mapping.

The original Xenopus canonical input has9,639proteins:8,649map to unique
reference numeric IDs and990do not. The mapping's XENTR interval contains
17,272IDs, of which8,623are absent from the original input. These are
identity-coverage counts, not independent gene-family observations.

## Scientific Consequences

This is broader than the14missing SwissTrees accessions. Original-input
inference cannot predict relations involving those absent reference IDs,
while predictions containing unrecognized accessions are removed by the
existing mapping filter. The audit does not quantify the effect on all
challenge scores or establish that all historical tools used identical
sequence inputs. It does not equate different accessions by sequence or
prove identity of scorer sequence bytes for recognized accessions.

Treat original-input QfO tables and current factorial results as explicitly
release-limited evidence until corrected-release compatibility is checked.
Do not fix this by renaming accessions, dropping inconvenient reference
relations, replacing active inputs or combining releases within one table.
Preserve original results. Corrected-input inference may change graph and
tree context beyond the replaced species, so a score-only substitution is
not an established remedy.

The separate corrected archive transfer21687 remains active. Its actual
canonical sequences, all78file differences and mapping coverage must be
verified before any new input freeze or rerun plan is admitted.

## Reproduction

```sh
python benchmark_tools/compare_qfo_corrected_archive.py \
  --archive qfo_benchmark/proteomes/QfO_release_2020_04.tar.gz \
  --prepared benchmark_tools/results/qfo_factorial_prepared_20260917.json \
  --mapping qfo_benchmark/benchmark-webservice/reference_data/2020/mapping.json.gz \
  --aliases benchmark_tools/results/swiss_sequence_alias_audit_20260917.json \
  --output benchmark_tools/results/qfo_original_archive_mapping_verified_20260917.json
```

The comparison script's name reflects its planned corrected-release use;
the report explicitly records the actual original archive identity. Its
generic status does not assert corrected-release provenance. The machine
report retains every canonical file hash, all990unmapped accession names,
all8,623missing numeric IDs and species counts. Sources are hash-checked
before/after reading; gzip is consumed through EOF. Thirteen tests cover
file inventories, identities, missing reference recovery and native species
interval boundaries. No archive files were extracted or overwritten.
