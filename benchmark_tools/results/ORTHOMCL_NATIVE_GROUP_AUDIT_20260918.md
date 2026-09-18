# Native OrthoMCL Final-Group Audit

## Scope

`audit_orthomcl_native_groups.py` independently checks `all_orthomcl.out`
against the native `tmp/all_ortho.mcl` partition, `tmp/all_ortho.idx` raw
identifier index and input `all.gg` ownership map. It does not use
pre-clustering graph edges as ortholog predictions.

The pinned OrthoMCL 1.4 `orthomcl_module.pm` implementation of
`mcl_backindex` (lines 776-830) omits clusters of size one and emits every
larger cluster, including those containing only one taxon. Final labels
retain the MCL cluster IDs; omission of singletons can leave gaps. The
auditor follows these semantics, without adding singleton groups.

Validation covers dense unique raw indexes, complete MCL partition coverage,
matrix dimensions and delimiters, wrapped cluster rows, unique membership,
exact final-cluster IDs and sizes, taxon annotations and counts, and omission
of singleton clusters only. All final groups of size at least two must be
present. An all-singleton partition may legitimately have an empty final
file. Raw gene identifiers are not converted to accessions at this stage.
Cross-species clique counts are calculated analytically, without expanding
all pairs. Input/source files are hashed before and checked after validation.

## Retained Evidence

All four retained native runs pass:

| Run | Input proteins | Indexed/grouped | Final groups | Single-species groups | Cross-species clique pairs |
| --- | ---: | ---: | ---: | ---: | ---: |
| Original serial fixture | 42 | 41 | 12 | 6 | 30 |
| Patched, one worker | 42 | 41 | 12 | 6 | 30 |
| Patched, two workers | 42 | 41 | 12 | 6 | 30 |
| Fresh staged inputs, two workers | 42 | 41 | 12 | 6 | 30 |

Each has one input protein absent from the native index, no singleton MCL
clusters, and one ungrouped input protein. These checks independently compare
each run's final groups to its own partition; the earlier three-arm probe
supplies the separate cross-run partition comparison.

Machine-readable reports and SHA-256:

- `orthomcl_native_group_audit_native_serial_20260918.json`:
  `a4bcefc25ee6524f8aba553c2b6036401b52e01bf402359dd617dc35a9032417`
- `orthomcl_native_group_audit_patched_one_20260918.json`:
  `f3e697d41e05a6695f1a5fbd0d1fb170392d8b44b6d667455b27bc65bb544aec`
- `orthomcl_native_group_audit_patched_two_20260918.json`:
  `f9b7038d80bba3b46829a10fc2d2bd2e0da59ea1e97c2bf9a55cdd97e91e9264`
- `orthomcl_native_group_audit_fixture_20260918.json`:
  `eda23e2e6ba13592599b3464646850311aed261542f7ddf19616564a4130af2a`

The four reports have status `native_final_groups_match_mcl_partition`;
`accuracy_admitted` and `publication_ready` remain false. They are component
audits, not scheduler or full-runtime admission records.

## Verification And Next Gate

51 new tests cover corrupt headers, dimensions, counts, IDs, taxa, duplicate
and missing members, singleton omission, raw IDs, wrapped rows, truncation,
trailing content, output provenance and refusal to overwrite outputs.
The focused suite, including the native input-staging test, passed 109 tests
in 8.69 seconds with `ORTHOHMM_LEGACY_BLAST_SMOKE=1`.

Production inference job 21750 remains pending behind corrected BLAST and
BPO preparation/admission. Do not modify its frozen executor. The next
independent admission must bind terminal scheduler evidence, execution
manifest, source/runtime identities, staged input preservation, complete
native outputs and this partition audit before conversion/scoring. Its
native final groups should be scored as cluster-derived cross-species
cliques, not native pairwise phylogenetic orthologs. Search failures and
ungrouped input coverage must remain visible. No corrected production result
has been admitted by this milestone.
