# Corrected BPO Preparation Submission

Submitted job **21748** with `afterok:21746`: it depends on successful
independent corrected BLAST admission, not merely a zero BLAST exit code.
`scontrol show job 21748 -o` confirms PENDING/Dependency, two CPUs,
64 GiB RAM, a 24-hour limit, node request bizon, no requeue and zero
restarts. No corrected BPO preparation has executed yet.

Frozen executor:
`benchmarks/work/publication_qfo_corrected_bpo_v1`, revision
`826e963cf06ed414609f0609b96c6bfd283fc2d8`.
Batch file: `benchmark_tools/results/qfo_corrected_bpo_batch_20260918.sh`
inside that executor. Native admission job argument:21746.

The batch verifies the clean executor revision, computes the completed
admission report hash and launches the dedicated interpreter with a clean
environment, isolated mode, disabled bytecode writes, and an absent,
job-specific bytecode-cache prefix. Numerical-library thread counts are one.
The wrapper checks the two-CPU/64-GiB allocation on bizon and the complete
search-admission provenance before preparing a fresh checkpoint.

## Runtime Enforcement

`verify_orthomcl_python_runtime.py` binds the retained manifest SHA-256
`3519dbf43a376c935cc76ab00559af045c3e94ce47bf421319b970e660c8cddc`,
requires the expected interpreter/environment/package identity and requires
every currently mapped file to match a pinned runtime file. It re-inventories
all2,928runtime entries. The preparation wrapper invokes this before and
after conversion/index validation. New or changed mapped libraries fail;
an unchanged subset is allowed because clean-locale startup can map fewer
locale files than the initial inventory.

Historical imported-source records in the runtime snapshot are not
misrepresented as the current executor identity: source revision and clean
worktree checks separately bind the revised helper code.

## Verification

63focused tests passed with legacy-runtime tests enabled. They cover Python
identity/package/library/manifest rejection, corrected admission binding,
and checkpoint success/failure paths. `bash -n` passed for the batch.

The frozen worktree passed a clean-environment runtime check, an actual
native checkpoint fixture and a subsequent runtime recheck. All six BPO
records, three query ranges and seven offsets were validated. `cmp` against
the native BioPerl fixture passed byte-for-byte. The executor remains clean.

Retained fixture report: `orthomcl_frozen_bpo_checkpoint_20260918.json`,
SHA-256 `96cb5bfc2259e2d5c350a9559bd9df25dc9c22598c2884a65ba34baca0444903`.
Raw artifacts: `benchmarks/work/orthomcl_frozen_bpo_checkpoint_20260918/`.

## Remaining Gates

The expected full-data checkpoint is under
`benchmarks/results/qfo_corrected_orthomcl_v1/bpo_preparation/`.
Its success status remains preparation pending independent admission.
Terminal scheduler verification and independent complete checkpoint
admission must precede native OrthoMCL inference. Final-group validation,
pair conversion, QfO scoring and failure-impact review remain outstanding.
Shared-host preparation accounting is not matched-resource timing evidence.
