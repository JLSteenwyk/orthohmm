# Corrected QfO Main Comparison

The [third corrected comparison export](qfo_corrected_comparison_20260918_v3/scores.md)
contains three admitted methods and five explicitly missing rows. Its JSON
manifest retains exact admission/conversion hashes, native precision and
recall, prediction semantics, pair accounting, exporter identity and helpers.
Earlier table versions and historical-input results are unchanged.

| Method | VGNC F1 | SwissTrees F1 | TreeFam-A F1 | Secondary mean |
| --- | ---: | ---: | ---: | ---: |
| OrthoHMM high sensitivity | 0.665933 | 0.685498 | 0.605008 | 0.689555 |
| SonicParanoid 2.0.9 | 0.982794 | 0.798459 | 0.771956 | 0.769506 |
| Proteinortho 6.3.6 | 0.954896 | 0.718111 | 0.643187 | 0.763216 |

The full export also includes GO similarity, EC similarity and FAS, which
are not F1. The mean is a project-defined secondary summary, not a primary
endpoint. These point estimates do not establish paired significance or a
complete corrected-release ranking.

## Publication-Method Identity

Only factorial p1_c0_r0 can supply the high-sensitivity row; only p1_c1_r1
can supply phylogeny satellite_v2. The latter is not yet admitted. The other
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
SonicParanoid contributes 15,248,739 native species-pair relations;
Proteinortho contributes 4,695,385 native post-clustering graph relations.
All have zero mapping losses. These counts describe prediction volume, not
protein coverage. Differing prediction semantics remain visible in both
human-readable and machine-readable exports.

The new adapter and both existing exporters pass 70 focused tests, including
rejection of nonpublication cells, duplicate methods, changed admission or
conversion hashes, unadmitted scores and invalid replay evidence. The actual
three-method export completed with the previously recorded admission hashes.
No parameter, score endpoint, frozen source helper or underlying prediction
was changed to assemble this table.
