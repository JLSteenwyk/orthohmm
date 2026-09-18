# Corrected FastOMA Input Staging

Preparation job **21738** is queued with `afterok:21731`, the independent
corrected OrthoFinder native-admission job. It reserves 2 CPUs and 64 GiB
for at most four hours on bizon. It does not launch FastOMA inference.

The immutable executor is
`benchmarks/work/publication_qfo_corrected_fastoma_prepare_v1` at
`82ef08e8a7d9511360ab52b62504871ba62753fc`.
The retained batch script is
`qfo_corrected_fastoma_prepare_batch_20260918.sh`.

The preparation checks successful admission-job accounting, the frozen
OrthoFinder admission source, the exact corrected primary plan, and the
admitted rooted node-labeled species tree. It requires 78 corrected FASTAs
and 984,137 proteins, including injective FastOMA accession transformation.
Tree leaves must exactly match the proteome species. The original-release
tree is not a fallback.

Fresh copies, not links, are written to
`benchmarks/work/qfo_corrected_fastoma_inputs_20260918/`:

- `proteome/*.fa`: byte-identical corrected source FASTAs.
- `species_tree.nwk`: byte-identical admitted corrected OrthoFinder tree.

The prospective manifest is
`benchmarks/work/qfo_corrected_fastoma_staging_20260918.json`.
It binds source and staged checksums, admission accounting, pinned FastOMA
assets and the completed tiny resource probe. Those records and copied
inputs are rechecked before success is recorded. Failed partial stages
are retained and never overwritten by a rerun.

Validation: 37 focused tests passed, including 15 staging tests; Bash syntax
validation passed. Actual production staging remains **pending**, because the
corrected OrthoFinder tree does not yet exist. Tests are not substituted for
that real execution. No inference or scoring result is implied.

Next steps: inspect the completed manifest, freeze the full fresh inference
command and runtime (including Java/Docker/Nextflow), then launch with the
prespecified resource configuration. Independently admit native outputs,
convert native pairs, and run corrected QfO scoring afterward. This comparator
uses a supplied OrthoFinder tree and must not be described as independent
FastOMA species-tree inference or dedicated matched timing.
