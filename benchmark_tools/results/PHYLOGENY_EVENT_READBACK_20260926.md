# Frozen Phylogeny Event and Pair Readback

Added an independent implementation of the frozen `species_overlap` root-group
rule and `positive_paralogy` pair rule in `derive_phylogeny_events.py`, with a
full-output reader in `audit_phylogeny_events.py`. Neither imports production
reconciliation code. Species mappings use DendroPy 5.0.8 MRCA queries; the
production implementation uses explicit ancestor-path intersections. Tests
include hand-checkable speciation, overlap, mapping conflict, branch support,
subroot duplication, supported satellites and unsupported satellite detachment.

The reader recomputes all node-table fields from saved rooted trees, including
parent identifiers, descendant genes/species, species mappings, event/pair-event
distinctions, confidence, overlap, mapping conflict and formatted support. It
checks reconciled-tree clades and event labels, derives pre-constraint root
groups and pair sets, then independently applies the recorded satellite
constraints. Final groups and complete native pair/confidence dictionaries
must match exactly. It does not assume that all root-group pairs are orthologs.

## Historical Full OrthoBench Result

All checks pass for retained p1_c1_r1:

| Quantity | Recomputed count |
|---|---:|
| Reconciled families | 8,681 |
| Bypassed families | 45,764 |
| Reconciliation node rows | 337,845 |
| Duplication events | 59,379 |
| Speciation events | 104,016 |
| Uncertain pair events | 13,034 |
| Final root groups | 59,770 |
| Final native ortholog pairs | 966,439 |
| Satellite constraints | 8,440 |
| Supported constraints | 5,852 |
| Detached constraints | 2,588 |
| Detached genes | 8,401 |
| Root groups added by membership filtering | 2,117 |
| Pairs removed by membership filtering | 16,885 |

The historical constraint sidecar resides in the factorial candidate directory,
not the final cell directory. Its 13,661,409 bytes match the previously frozen
preparation manifest: SHA-256
`55dee278b7feac6cc480c4e88b21adabc7ba093d8c9534501794e7efcb3270a9`.
The [tracked summary](installed_ob_historical_events_20260926.json) links this
provenance to the full local report, 34,745 checked records, SHA-256
`8bcc41ab29b5e84ea4fdbfd8e772c73f8a904bf3900a3a48c0ffb065422161c3`.
The large generated report remains local rather than being committed.

The installed synthetic fixture also passes;
[receipt](installed_phylogeny_events_fixture_20260926.json). All 85 related tests
pass. Integration mutations that pass structural validation but fail this
oracle include an omitted pair with matching summary counts, changed confidence,
an incorrect group split, changed node semantics and a hash-consistent incorrect
annotated-tree event.

```bash
python -m benchmark_tools.audit_phylogeny_events \
  --directory benchmarks/results/publication_ob_factorial_v1/cells/p1_c1_r1/orthohmm_phylogeny \
  --structure benchmarks/work/installed_ob_historical_structure_20260926.json \
  --constraints benchmarks/results/publication_ob_factorial_v1/candidates/p1_c1/orthohmm_working_res/phylogeny_candidate_merges.json \
  --output benchmarks/work/installed_ob_historical_events_new.json
```

## Interpretation and Remaining Scope

These checks establish consistency with the frozen rules conditional on the
saved rooted trees and recorded constraint sidecar. They do not establish
biological truth, optimal alignments/rooting, search correctness, candidate
selection or whether the reconciliation/bypass selection was appropriate.
The reader shares the DendroPy parser family with production. Hierarchical
group-table validation remains separate. Existing structural and residue
readbacks provide complementary checks, not independent biological evidence.

No scores, scientific source or running executor changed. Job 22179 still
requires successful native completion, frozen-provenance/input checks and all
readback stages before a fresh reproduction claim. Publication readiness and
overall superiority are not claimed.
