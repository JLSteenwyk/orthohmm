# Explicit High-CPM Candidate Handoff

This prospective execution amendment follows the separately admitted recovered
seed and its independent readback. It is not an accuracy-driven parameter change,
new optimizer attempt, or claim that the historical failures have been repaired.
The original admission 22155 remains failed and old candidate 22156 remains
cancelled. The default candidate entrypoint retains the original admission gate.

## Fixed Inputs And Science

Explicitly select `qfo_cpm_helper_recovery_admission_readback_20260930.json`, SHA256
`5fb31ea6cee69501d5ad1826789e9a178f230bb909012537e123d18d21d65850`.
Its admitted seed report is SHA256
`35febb4c1810892988beb79bc8e6ac2c6a2e1ffbf1b111492ad1649f39e04997`.
The handoff verifies both, the five source Git bindings, full transitive input
inventory, four stage memberships, and current scheduler success/failure identities.
The earlier admission streamed the entire graph; its independent readback used
a separate NumPy reader. Do not rerun native refinement or graph optimization.

Use the original refined seed from `qfo_cpm_checkpoint_recovery_v1`, not a new
seed selected from scores. CPM resolution remains 0.12, seed 4; final seed has
984,137 genes in 390,845 groups and SHA256
`f4c6f1973bc9636828baf1e6d9be3f416a18fc8ada5302fb502c082495fc1811`.
Preserve original missing profile counters and successful stage timings.

Candidate construction reuses the unchanged `prepare_qfo_cpm_candidates.py`
builder, SHA256
`b6f3325e2e6c33eecbf74256e2bd88d339015188e47c9296fedea1155029b978`.
Use `controlled_expansion` with label `control`, the frozen `satellite_v2`
parameters from `qfo_corrected_factorial_v1/manifest.json`, and the original
numeric hit checkpoint. No threshold sweep, new search, default promotion,
reference-label access or endpoint change occurs. Preserve the prespecified
seven-arm parameter family and its downstream multiplicity treatment.

## Runtime And Resources

The helper-complete runtime amendment authorized independent refinement only.
It does not authorize switching candidate science to that runtime. Candidate
construction must use `/home/bizon/anaconda3/bin/python`, the original frozen
QfO launcher/core and native libraries, and pass the unchanged runtime verifier
before and after construction. Its report must equal the corrected plan and
both baseline runtime records. If this comparison fails, preserve failure and
stop; never weaken it or silently substitute the new environment.

Use one new detached executor at the committed source revision on local bizon:
2 CPUs, 64 GiB, 4-hour Slurm limit, one task, no GPU, no requeue or automatic
retry. Use PYTHONHASHSEED=0 and OMP/OPENBLAS/MKL thread counts of one, disable
user-site packages and bytecode writes. Preserve the submission script,
command, source identities, scheduler state and stdout/stderr. Do not release
the cancelled old dependency or modify any unrelated job/service.

## Output And Admission

Only the all-three-options explicit CLI path may consume the new seed readback:
`--helper-recovery-readback`, `--helper-recovery-readback-sha256` and
`--helper-candidate-protocol-sha256`. Partial selection must fail before any
scientific imports or output writes. The protocol's SHA is supplied explicitly
and recorded. The default historical path must continue rejecting failed 22155.

Write a fresh `benchmarks/results/qfo_cpm_helper_recovered_candidates_v1` tree;
refuse existing outputs. Retain the seed and full admission lineage, frozen
settings/runtime/input identities and candidate content audit. Recheck inputs,
numeric checkpoint and runtime after construction. Failures remain failures;
do not silently resume, overwrite or retry them.

A successful build is `recovered_cpm_candidates_prepared_pending_admission`,
with explicit `seed_handoff=explicit_helper_runtime_seed_amendment`. It is not
candidate admission or accuracy. Independent candidate admission, phylogenetic
inference/reconciliation, pair conversion, scoring and robustness integration
remain separate required gates. No downstream jobs may bypass those gates.

This shared-host incremental construction does not measure comparable full
pipeline time or controlled efficiency. Controlled timing remains deferred.
