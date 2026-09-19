# Corrected QfO Main Comparison

The [fourth corrected comparison export](qfo_corrected_comparison_20260919_v4/scores.md)
contains four admitted methods and four explicitly missing rows. Its JSON
manifest retains exact admission/conversion hashes, native precision and
recall, prediction semantics, pair accounting, exporter identity and helpers.
Earlier table versions and historical-input results are unchanged.

| Method | VGNC F1 | SwissTrees F1 | TreeFam-A F1 | Secondary mean |
| --- | ---: | ---: | ---: | ---: |
| OrthoHMM high sensitivity | 0.665933 | 0.685498 | 0.605008 | 0.689555 |
| OrthoHMM phylogeny satellite_v2 | 0.901690 | 0.833513 | 0.614864 | 0.761510 |
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
