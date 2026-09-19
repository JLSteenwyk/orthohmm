# Corrected QfO Main Comparison

The [fifth corrected comparison export](qfo_corrected_comparison_20260919_v5/scores.md)
contains six admitted methods and two explicitly missing rows. Its JSON
manifest retains exact admission/conversion hashes, native precision and
recall, prediction semantics, pair accounting, exporter identity and helpers.
Earlier table versions and historical-input results are unchanged.

| Method | VGNC F1 | SwissTrees F1 | TreeFam-A F1 | Secondary mean |
| --- | ---: | ---: | ---: | ---: |
| OrthoHMM high sensitivity | 0.665933 | 0.685498 | 0.605008 | 0.689555 |
| OrthoHMM phylogeny satellite_v2 | 0.901690 | 0.833513 | 0.614864 | 0.761510 |
| OrthoFinder 3.1.5 full | 0.988546 | 0.848413 | 0.791918 | 0.787663 |
| OrthoFinder 3.1.5 sequence-only checkpoint | 0.142755 | 0.690515 | 0.734407 | 0.548352 |
| SonicParanoid 2.0.9 | 0.982794 | 0.798459 | 0.771956 | 0.769506 |
| Proteinortho 6.3.6 | 0.954896 | 0.718111 | 0.643187 | 0.763216 |

The full export also includes GO similarity, EC similarity and FAS, which
are not F1. The mean is a project-defined secondary summary, not a primary
endpoint. These point estimates do not establish paired significance or a
complete corrected-release ranking.

## Publication-Method Identity

Only factorial p1_c0_r0 can supply the high-sensitivity row; only p1_c1_r1
can supply phylogeny satellite_v2. Both are now admitted. The other
six factorial cells are ablations and are rejected as publication-method
substitutes. In particular, the p0_c0_r0 initial-HMM sequence-control score
must not be relabeled as the full high-sensitivity method.

The exporter requires the pinned [independent replay admission](QFO_CORRECTED_REPLAY_ADMITTED_20260918.md)
that established exact corrected native/replayed high-sensitivity partition
equality across 984,137 genes. It then reuses the existing comparator and
factorial admission adapters for status, identity, conversion, prediction
semantics and native endpoint arithmetic. Frozen scoring/uncertainty helpers
are unchanged. The export does not repeat upstream inference or its audits.

High sensitivity contributes 9,032,719 cross-species group-derived pairs;
phylogeny satellite_v2 contributes 5,959,560 native phylogenetic pairs;
SonicParanoid contributes 15,248,739 native species-pair relations;
Proteinortho contributes 4,695,385 native post-clustering graph relations.
All have zero mapping losses. These counts describe prediction volume, not
protein coverage. Differing prediction semantics remain visible in both
human-readable and machine-readable exports.

The new adapter and both existing exporters pass 70 focused tests, including
rejection of nonpublication cells, duplicate methods, changed admission or
conversion hashes, unadmitted scores and invalid replay evidence. The actual
four-method export completed with the previously recorded admission hashes.
No parameter, score endpoint, frozen source helper or underlying prediction
was changed to assemble this table.

The fourth export adds final scoring admission 21788 with SHA-256
`49b7d837b2ba2928b0676c9974e2ec16db5211f1086c366c7bc11805a2ee751f`.
Its manifest SHA-256 is
`b8f37d88b16e64b34a8ec01803734459a7dadef59272dec985d0550c1735ca67`.
All three previously admitted rows compare equal as complete structured
records, and all recorded input/helper/output hashes were rechecked.
SwissTrees phylogenetic precision/recall are 0.955177/0.739341; TreeFam-A
precision/recall are 0.959110/0.452464. These expose different precision-recall
trade-offs; greater SwissTrees F1 than another admitted method is only a
point-estimate comparison, not a tested superiority claim.

## OrthoFinder Admission And Fifth Export

Both corrected OrthoFinder score validators completed successfully:
`21736_0` in 3:02 and `21736_1` in 3:28 (two CPUs on `bizon`, exit 0:0).
These are validation durations, not inference timings. The full reports
remain under `benchmarks/work/` (approximately 146 MiB each); they are not
committed as large raw artifacts.

| Admission report | SHA-256 |
| --- | --- |
| `qfo_corrected_orthofinder_full_assessment_admission_20260918.json` | `f0c7064413622db6ebe02965d8781b5c6862b5246ad10cdd6b907b3c1815a46c` |
| `qfo_corrected_orthofinder_sequence_only_assessment_admission_20260918.json` | `ecd8721dc2e8f2399ccf86fd7d89f917447325ee4e812bae29bfaf3b2bbe851e` |

Full OrthoFinder contributes 14,215,382 native phylogenetic pairs;
sequence-only contributes 163,277,439 cross-species clique pairs from its
pre-phylogenetic MCL groups. Both have zero mapping losses. The latter is
a diagnostic checkpoint, not native pairwise orthology output, and its
pair volume must not be called protein coverage.

Full OrthoFinder has higher point-estimate F1 than phylogenetic OrthoHMM
on VGNC, SwissTrees and TreeFam-A. OrthoHMM has higher GO/EC similarity
and FAS. SwissTrees precision/recall are 0.937853/0.774548 for full
OrthoFinder versus 0.955177/0.739341 for phylogenetic OrthoHMM. TreeFam-A
precision/recall are 0.924527/0.692578 versus 0.959110/0.452464. These
observations describe a precision-recall trade-off; they do not establish
its mechanism or paired statistical significance.

The sequence-only checkpoint has VGNC precision/recall 0.076891/0.995362,
SwissTrees 0.569068/0.877863 and TreeFam-A 0.620470/0.899600. Reporting only
its high recall or its secondary mean would obscure that trade-off.

The unchanged exporter completed using all six pinned admission hashes.
All four previous method rows are identical as complete structured
records, and all export input/helper/output hashes were independently
rechecked. All 70 exporter tests pass. The fifth manifest is 19,317 bytes,
SHA-256 `7256920214d2eeda9cf596c7c39212992c7f3fcc36dd61e1848fa44d5725e01e`.
Its `checked_records` and per-method `admission` fields give exact inputs
for reproducing the export with `export_qfo_complete_comparison.py`.

FastOMA and OrthoMCL remain unadmitted in this corrected-release table.
Paired comparator uncertainty and controlled resource comparisons remain
unfinished; historical intervals cannot be attached to these corrected
point estimates. No parameters or endpoints changed after viewing scores.
