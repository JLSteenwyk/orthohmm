# Corrected Proteinortho Command Freeze

Prepared, not submitted. This preserves the original invocation from
`qfo_benchmark/run_tool.slurm`: Proteinortho 6.3.6, project `qfo`, 32 CPUs,
all 78 canonical FASTAs and otherwise default scientific settings.
Native graph semantics remain the post-clustering `qfo.proteinortho-graph`,
not the pre-clustering graph and not clique-expanded orthogroups.

Manifest: `qfo_corrected_proteinortho_commands_20260918.json`, SHA-256
`65f9baedeb6f6cbb66198e2a6a24ae8531b48eea9a04809c2c2443d570d4f93a`.
Image: 131,035,136 bytes, SHA-256
`990be9d066a02302fd49520abfa28d63bfdb5ad0c4792a99ac75fc6152cbe561`.
Read-only probes all exited zero: Proteinortho 6.3.6, DIAMOND 2.1.12,
Singularity-CE 4.3.2. Proteinortho and DIAMOND versions agree with the
retained original inference log. Current resolution does not prove
historical executable identity.

The preparer checked the pinned comparison registry, corrected primary
manifest and all corrected inputs, ran the existing primary environment
verification and bound original launcher/converter sources. It does not
copy inputs, run inference, convert graphs or evaluate scores. Its output
is explicitly unauthorized for execution pending launcher preparation.

Validation: nine command-construction tests passed, including wrong
counts, duplicate names, path traversal, option-like names and unsupported
extensions; 27 tests passed when combined with the corrected primary
planner/runner tests. The real preparation completed successfully and
rechecked recorded files after version probes.

Next: pin the container runtime configuration, inherited environment and
bind-mounted dependencies; implement fresh-copy execution with pre/post
checks and failure recording; freeze conversion commands; then submit
under shared-host descriptive resource accounting. Do not run the older
wrapper that deletes its staging directory. Other corrected comparator
rows, factorial reruns and scoring remain required.
