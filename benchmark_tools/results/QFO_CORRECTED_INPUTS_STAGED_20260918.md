# Corrected QfO Inputs Staged And Inventoried

Native sequence audit21689completed0:0 in3:14 with pinned executor
`25394c3e2b5b7edf92b710d5edb6859b1d09866c`. Report
`qfo_corrected_archive_sequences_20260917.json` has SHA-256
`b939e03bfaa4975e0a3b38bbedfcd39f9ab18eff05e956cdc6967dce23c065e1`.
It covers all984,137native identities:983,959exact sequence matches and178
differences fully explained by B/O/U/Z-to-X representation. No unexplained
difference or absent native identity remains. The178representation differences
span24species, including5XENTRentries; exact IDs and hashes remain in the
report. Representation equivalence is not exact byte identity.

After reviewing both successful scheduler audits, their pinned clean
executors, common archive identity and empirical results, ran the stager
from detached executor `publication_qfo_corrected_staging_v1` at
`1096ccea2cdd5f7cbbab7dec5ed69282cca00166`. Supplied both reviewed report
hashes explicitly. Original inputs and original-release runs were untouched.

New canonical directory:
`benchmarks/work/qfo_corrected_inputs_20260918/`.
The allowlisted78FASTA files were extracted unchanged; no DNA/additional
files, isoform substitutions or residue normalization. Staging manifest
SHA-256 is
`07a890eb816944f946d46039559a6046c2a3664f033eaa0f30f35bb25b9c9ab8`,
preserved as `qfo_corrected_staging_manifest_20260918.json`.

The same frozen executor then ran `audit_qfo_corrected_staging.py`, binding
that exact staging-manifest hash and the retained reference mapping. Direct
inventory confirms78distinct reference species/proteomes,984,137unique
sequences and440,246,934residues. Every numeric species interval is completely
covered, with no mixed-species proteome, duplicate ID or unexpected file.
All staged hashes were checked before and after parsing.
Report `qfo_corrected_staged_inventory_20260918.json` has SHA-256
`c86c2d6337928938a5a682de0761cc0d293143a1fc5760fe39e5e0c8065e97ca`.

Both reports retain `inference_authorized: false`. Next, prepare and commit
the immutable execution manifest covering exact tool/runtime identities,
original scientific settings, inputs, commands, resources, seeds, conversions
and separate output namespaces. The corrected-release protocol continues
to require all eight comparator rows and the eight factorial cells. Do not
reuse original groups, inferred trees or scores as corrected-input outputs.
This staging milestone establishes input compatibility, not corrected
accuracy, independent biological validation or publication readiness.

Fifty-eight targeted staging, direct-inventory and assessment-admission
tests pass at this milestone. Production staging and direct inventory both
completed successfully; neither launched inference. Large FASTAs remain
outside git, with committed content/provenance manifests.
