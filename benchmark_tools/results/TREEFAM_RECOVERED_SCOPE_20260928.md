# Recovered Subtrees Cannot Cover The Complete Reference

The recovered Selectome archives contain only Euteleostomi subtrees. To test
whether that scope could cover the retained QfO TreeFam reference, this audit
uses the existing QfO identifier mapping and retained lineage metadata, without
mapping any Selectome gene or inventing original-family labels.

## Necessary Coverage Exclusion

| Reference quantity | Count |
| --- | ---: |
| Proteins incident to retained relations | 11,130 |
| Incident proteins outside Euteleostomi | 5,334 |
| All retained relations | 79,320 |
| Relations with at least one outside-clade endpoint | 55,933 (70.52%) |
| Outside-clade-endpoint ortholog relations | 35,868 |
| Outside-clade-endpoint paralog relations | 20,065 |
| Relations with both endpoints inside Euteleostomi | 23,387 |

Outside-clade proteins comprise 1,241 Arabidopsis, 1,106 Caenorhabditis,
706 Ciona, 1,184 Drosophila, 502 Schizosaccharomyces and 595 Saccharomyces
records. The [complete audit](selectome_reference_scope_20260928.json) includes
all reference species counts, full scientific names and input hashes.

Thus at least 70.52% of the retained reference relations cannot be represented
by the recovered vertebrate-only subtree collection under a species-consistent
mapping. This is a necessary exclusion based on taxonomic scope, not a measured
reconstruction rate or an accuracy-effect estimate. The 23,387 inside-clade
relations are only potentially in scope; gene, family and event coverage for
them remains unestablished. No subgroup score or confidence interval was made.

## Evidence And Checks

The prior native count audit establishes 79,320 relations and 11,130 incident
proteins; its retained raw export is checked and parsed again. QfO's existing
mapping translates those raw reference identifiers to its own entry numbers.
The native offset convention identifies their species. Entry assignments are
complete and injective for the incident reference proteins. The retained
phyloXML defines one Euteleostomi clade; every observed reference species has
taxonomy metadata in that file.

The phyloXML title names QfO 2018, despite residing in the retained 2020 resource
directory. It is used here as taxonomy metadata, not assumed to be a new
2020 species-tree inference. Ten reference members with no retained relation
are excluded from this incident-protein analysis and explicitly counted as
unanalyzed. Ortholog/paralog classes are reconstructed from TP/FN versus FP/TN
labels; these counts do not compare prediction quality.

All 17 focused scope/count-reader tests pass, including counting a relation
with two outside endpoints only once and rejecting missing, extra or empty
species assignments. Source and retained evidence hashes are checked, and the
audit fails on inconsistent reference counts or species coverage.

```bash
python -B -m benchmark_tools.audit_selectome_reference_scope --repo . --output /tmp/selectome_reference_scope.json
```

The existing QfO mapping is not `treefam2reference.txt`: it is used only to
identify species of already retained reference proteins. The original mapping
and complete source trees remain missing. Recovering vertebrate subtrees does
not resolve full-reference family uncertainty; no benchmark inputs, scores or
defaults changed, and no external contact or DGX work occurred.
