# Conditional Corrected-Release Evaluation Protocol

## Status And Purpose

Written before the corrected archive has finished downloading and before any
corrected-input predictions or scores exist. This fixes the evaluation policy,
not an empirical input admission or executable inference manifest. Acquisition
21687 and dependent audits 21688/21689 remain the prerequisite workflow.
No extraction, inference, scoring or overwrite is authorized by this document.
QfO is development-exposed; correcting its input release does not turn it into
independent confirmation or permit retuning the method.

Original retained comparator inputs agree at the evidence levels recorded in
`QFO_COMPARATOR_INPUT_PARITY_20260917.md`. The original archive is nonetheless
incompatible with the scorer's complete Xenopus reference. Uniform input
limitations are not evidence of uniform effects on accuracy or rankings.

## Input Acceptance Before Any Rerun

1. Require successful terminal acquisition and both audit jobs, inspect their
   logs, and verify their pinned executor revisions and source identities.
   Scheduler success alone is insufficient. Recheck the actual archive length
   2,648,666,198 bytes and SHA-256 against both audit input records. Preserve
   source URL, HTTP metadata, logs and local checksum; do not describe the
   latter as a publisher-provided checksum.
2. Require agreement between both audits on the complete 78-file canonical
   inventory, archive member names, file sizes/hashes and all mapping counts.
   Inspect gzip EOF validation. The expected correction is confined to
   `UP000008143_8364.fasta`; any other changed canonical file requires a
   documented release-scope review before acceptance, not silent approval.
3. Require complete unique coverage of all 984,137 scorer numeric identities,
   no unmapped canonical input accessions, and recovery of all 14 previously
   absent SwissTrees accessions at their exact reference numeric identities.
   Do not drop truth relations, substitute isoforms or infer missing aliases.
4. Require complete native database coverage and no unexplained sequence
   differences. Report exact matches separately from B/O/U/Z-to-X-only
   representation matches; inspect their IDs and species. Preserve original
   archive residues for every method. Do not normalize inputs to make this
   gate pass or count representation matches as exact byte matches.
5. Freeze a separately named canonical input directory and manifest only after
   those empirical checks. Extraction must use the explicit audited regular
   member allowlist into a new directory, verify every output size/hash, and
   exclude DNA/additional files. Never run the destructive original preparation
   script on the existing input directory. Record counts, ownership and
   sequence/ID uniqueness again after staging.

If any condition fails, retain the reports and investigate the discrepancy.
Do not relax a requirement based on a favorable tool score. Record any needed
protocol amendment before predictions or scores are inspected. Compatibility
with retained scorer sequences does not independently validate annotations.

## Frozen Scientific Comparisons

Retain all eight comparator rows: OrthoHMM high-sensitivity and satellite_v2;
OrthoFinder 3.1.5 full pipeline and sequence-only MCL checkpoint;
SonicParanoid 2.0.9; Proteinortho 6.3.6; FastOMA 0.3.5; OrthoMCL 1.4.
Use the recorded original scientific settings, frozen OrthoHMM core 7f3a9e4
and source-equivalent launchers. Do not promote the later development banding
change. Freeze exact executable/container hashes, arguments, environment,
seed, resources, output semantics and conversion commands before launch.
Version labels alone are insufficient provenance.

The comparator registry is
`publication_comparison_orthomcl_complete_20260916.json`, SHA-256
`094f842ad8211450519d4238c0b3a5020a7969049d4b7306f5465d5acbf914a3`.
Satellite configuration record:
`phylogeny_satellite_v2_method_freeze_20260902.json`, SHA-256
`90f3f496211451bcc2ed4bd5a7d51736e8755dd719e1e10c208ed2e157941802`.
These are configuration/provenance anchors, not a claim that a recovered
OrthoHMM replay reproduces every historical partition. Any such difference
must remain explicit; no score transfer between nonidentical partitions.

Derive the sequence-only OrthoFinder row from the corrected full run's audited
MCL checkpoint using the same conversion as the original diagnostic. Preserve
native pair predictions for full phylogenetic QfO rows rather than substituting
RootHOG cliques. FastOMA retains the supplied-tree diagnostic design but must
receive the newly inferred corrected-input OrthoFinder tree, not the old tree.
It is still not an independent end-to-end tree inference comparison.

Rerun the eight P/C/R factorial cells from corrected-input HMM evidence with
the existing design in `QFO_FACTORIAL_PROTOCOL_20260917.md`, SHA-256
`f8946e12cefcf84abbee0fb9492f240c05508e045efe00a3006304d34c1fd115`.
P-off still contains initial HMM search. R-on scores native inferred pairs;
RootHOGs serve integrity checks only. Preserve original factorial jobs/results
as the original-release experiment; do not interrupt them or reuse their
scores as corrected-release measurements.

## Reuse And Failure Rules

Changed proteome content can affect candidate ranking, normalization, groups,
profiles and trees throughout the dataset. Recompute all affected downstream
inference; merely rescoring old predictions or appending recovered proteins
is not a valid corrected-input rerun.

Low-level search work for unchanged sequence pairs may be reused only after
exact input, program, parameter and output identity checks establish that its
values do not depend on changed database size/composition or other global
state. Recompute statistics when their database dependence requires it. In the
absence of that proof, rerun the search. Keep reuse records explicit, and
never label reused-stage wall time as end-to-end runtime.

Record all failures, including sequence-specific BLAST errors. Preserve the
OrthoMCL legacy engine/settings; isolate failed queries and quantify their
reference impact rather than silently dropping them or changing engines.
Correctness-driven reruns on bizon are not matched timing evidence. Do not
compete with the ongoing dedicated DGX timing array. A changed input release
does not relabel the existing scaling dataset or timing measurements.

## Scoring And Reporting

Keep the scorer reference, containers and commands pinned by
`qfo_assessment_environment_20260917.json`, SHA-256
`e86545fd04cb644ed642cf2a31b4993fb2225ca4c05743cd03e661db129e82bc`.
Use separate corrected-release prediction IDs, work directories, outputs and
admission reports. Verify native completion, conversion semantics, unique
cross-species pairs, mapping losses and benchmark task success per row.

Report all six QfO endpoints, available precision/recall and prediction
coverage. Keep the six-metric mean explicitly secondary. Do not remove a
method, family or metric because its corrected result is unfavorable.
Original and corrected results must be separate panels with clear input
checksums. Paired release differences are descriptive unless their own
prespecified dependent-data inference is added before score inspection.

For corrected SwissTrees comparisons retain the existing 18-family paired
bootstrap protocol, 100,000 PCG64 draws, seed 20260920, eight contrasts and
24-endpoint Bonferroni family. Protocol SHA-256:
`a1ba79268b732824e0df9c56dc1b14dce7b9f33c955523e8b87133c5440161b8`.
Reconstruct corrected raw confusion counts and recompute native raw/2+1
macro precision/recall and harmonic F1 within every draw. Do not carry over
original counts or intervals. The corrected factorial retains its separate
42-endpoint uncertainty family and seed 20260922. These are explicitly
separate analyses, not global multiplicity control across all publication
results. Do not interpret a change in statistical significance as a
significant release-by-method interaction.

This protocol is incomplete operationally until actual input admission,
executable manifests and resource scheduling are committed. It does not
declare corrected results available or the publication goal complete.
