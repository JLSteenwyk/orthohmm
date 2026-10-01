# Explicit Recovered High-CPM Inferred-Phylogeny Handoff

Prospective execution amendment, October 1, 2026. No new high-CPM phylogeny
or accuracy has been evaluated under this amendment. This is a separate
recovery-aware entry point, not relaxation of the historical failed gates.

## Admitted Inputs

Construction 22385 and independent candidate admission 22386 are COMPLETED
0:0 on bizon, each 2 CPUs/64 GiB. The full candidate admission is
`benchmarks/work/qfo_cpm_helper_candidate_admission_20260930/status.json`,
SHA256 `55c696a037d13b48d47dc752fcb6b69b5b8bba6a7b990acfb08b0d8fbf6e07c0`.
Its independent readback is
`qfo_cpm_helper_candidate_admission_readback_20261001.json`, SHA256
`8a5e832febbccf7d4d55a5d2b16e4454a54f83fc3c5570d9c16dfed6b7f12e3a`.
All 11,513 bound records and complete memberships/merge reconstruction were
validated; no accuracy was evaluated. Retain the full inventories locally.

The c3f1968b913a83ba47e7fb021347bdd83e2f78fc admission executor and source
SHA256 `cc9745c1f1cb8ddf4c595c6cc85b1cd7fec8fce1d68e8c1bd6e4cd657f68786f`
must remain unchanged. Before native inference, run that same frozen read-only
candidate admission in a separate fresh process and require exact JSON equality
with the retained report. This does not reconstruct candidates again.

Use only cpm_high/index 1: 984,137 genes, 390,845 seed families,
346,866 candidate families and 43,979 merges. Candidate partition SHA256
`def5d6f8d77b043745845774cabf2832757ed61951b87aace0ef102a0950613c`;
membership constraints SHA256
`894881980ba25c3526197fdb77464fa602f859d8f48f8cecd11486e8877b9d0b`;
recovered seed SHA256
`f4c6f1973bc9636828baf1e6d9be3f416a18fc8ada5302fb502c082495fc1811`.
Retain all four candidate sidecars and their original exact expansion parameters,
satellite_v2/high_confidence_pair policy and explicit helper-runtime seed amendment.

## Frozen Inference

Reuse `run_qfo_parameter_phylogeny.verify_baseline` unchanged. It admits the
corrected p1_c1_r1 baseline, original candidate/native outputs, full scientific
environment and package/tool/source inventories. A mismatch is a failure, not
authorization to rebuild an environment or substitute the helper runtime.
The helper runtime amendment applied to independent seed refinement only.

Use the unchanged `variant_cell` and `native_command` helpers. Only the admitted
candidate/constraint files and new destination replace baseline arguments;
append the original validated baseline as `--checkpoint-source`. Preserve
all full-pipeline scientific settings: 32 CPUs, MAFFT/FastTree, reconciliation,
root/pair rules and species-tree rooting. Require `--species-tree-mode infer`
with no supplied species tree or reference/competitor tree, and no benchmark
scoring arguments during inference.

This is incremental inferred phylogeny with validated checkpoint reuse, not a
new end-to-end sequence search or all-tree run. The frozen pipeline checks
matching family content/configuration and raw-tree hashes before reuse; eligible
species-tree inference checkpoints are also governed by its existing checks.
Retain actual checkpoint-hit counts in outputs. Do not call any of these
measurements controlled timing or pool them with the deferred scaling panel.

## Execution And Failure Policy

One new standalone local bizon task: 32 CPUs, 192 GiB, 24 hours, no GPU,
no automatic requeue or retry. This matches the historical native CPM
phylogeny allocation, not the separate prospective controlled timing budget.
Fresh destination: `benchmarks/results/qfo_cpm_helper_phylogeny_v1/cpm_high`.
Refuse existing directories or symlinks; never overwrite a prior attempt.
No dependency release, resubmission or reinterpretation of cancelled 22156
or failed original seed admissions 22155/22158/22082_1.

Use `/home/bizon/anaconda3/bin/python` and the original full-pipeline environment,
including its frozen user-site inventory. Pin hash seed 0 and OMP/BLAS/MKL
thread controls to 1. A clean detached driver executor must be committed and
pushed before submission. Verify the supplied protocol SHA256 and all source/
runtime/input identities before inference and again after native success.
Bind actual argv, launcher source equivalence, resolved tools, cwd, source,
job and file inventories. Preserve stdout/stderr, time log, execution status
and error postflight; restore cwd even on failure. No unrelated job, service,
DGX or workload-isolation change is authorized or needed for this accuracy arm.

## Admission And Interpretation

Successful execution is only `complete_pending_native_validation` with
accuracy_evaluated/native_outputs_validated/publication_ready all false.
Native output validation, lossless pair conversion, score admission and the
prespecified seven-arm/18-endpoint parameter analysis remain separate work.
No parameter, endpoint, multiplicity or publication-default change is permitted
based on this recovery. Original missing optimizer statistics and runtime
failure cause remain unresolved; do not infer accuracy or general superiority
from validated candidate memberships or completed phylogeny alone.
