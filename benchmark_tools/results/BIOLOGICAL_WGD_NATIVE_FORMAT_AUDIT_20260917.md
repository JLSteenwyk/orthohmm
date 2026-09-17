# Native Application Output Audit

All four tasks21661_0..3completed0:0. Native output integrity is admitted;
biological endpoint comparisons have not yet been calculated. This note
clarifies extraction mechanics without changing the frozen root-HOG endpoint.

## OrthoFinder 3.1.5 Root Table

The initial comparator validator expected the public
`Phylogenetic_Hierarchical_Orthogroups/N0.tsv` and rejected its absence.
This was an extraction assumption failure, not a failed inference run.
Installed source explains the lifecycle:

- `gene_tree_inference/trees2ologs_of.py`, lines320-335, constructs parallel
  original-ID and internal-ID membership rows, emitting the latter asN0.ids
  whenever the former is emitted for rootN0. SourceSHA256
  `e1e570832d85d807d16579718648b735e792909c07a5a69012125b49ea600d50`.
- `comparative_genomics/orthologues.py`, lines679-683, removes publicN0.tsv
  and legacy outputs whenrm_legacyis enabled. SourceSHA256
  `64a4a11a0a184572c14bfbc71eef977f5cd7ade7ed41adaf1f85ae93e3354f65`.
- `file_updates/ogs.py`, lines31-38, constructs finalOGs from rootHOGs and
  adds unassigned genes as singleton sets. SourceSHA256
  `e45c16079aa1188f0784323075eaeac09e8fd3c677d8e4884dc83f58720e3911`.

Use retained `WorkingDirectory/N0.ids.tsv` with exact `SequenceIDs.txt`
restoration. Require a bijective mapping over the complete input universe,
correct species columns, unique memberships and the four-species rootN0 tree.
The restored non-singleton groups must equal final native non-singletonOGs.
Do not supplement missing root memberships with the final-output singletons:
the prespecified endpoint remains explicit rootHOG assignment, not synthetic
singleton separation. No manual gene remapping or rerun was performed.

Actual root table:5,581groups,23,233assigned proteins,637unassigned proteins.
The5,925-group MCL checkpoint remains a separate diagnostic. Final admission
is `biological_wgd_orthofinder_admission_v2_20260917.json`; the earlier successful
report is retained but superseded after an unrelated Sonic log-parser fix
changed the shared validator source hash.

## SonicParanoid Log and Groups

Initial admission rejected the repeated identical `Main output directory`
line printed by two native stages. The parser now permits consistent path
repetition but rejects conflicting paths, repeated run headers, changed
settings or absent completion. A synthetic regression test covers both cases.
The fixed2.0.9/default/32-thread configuration, DIAMOND very-sensitive mode and
non-graph-only setting are confirmed. No method settings were changed.

Actual native table:5,467groups,22,347assigned proteins,1,523unassigned proteins.
Species filename columns are checked explicitly; missing assignments remain
missing. Admission: `biological_wgd_sonic_admission_20260917.json`.

## Interpretation

Counts above establish output integrity and coverage, not accuracy or biological
superiority. Preserve all240experimental pairs, the fixed eligibility counts,
and all six prospective examples in the forthcoming score report. Native
cross-species pair tables are not substituted for same-species anchor grouping.
