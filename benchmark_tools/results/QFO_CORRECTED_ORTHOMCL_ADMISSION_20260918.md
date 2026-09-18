# Corrected Native OrthoMCL Admission

## Gate

`admit_qfo_corrected_orthomcl.py` validates completed native inference from
the frozen runner `22f4eec4314304a0cda71e02595adfb20827c34b` at
`benchmarks/work/publication_qfo_corrected_orthomcl_v1`. It does not alter
that executor or restart inference. It requires:

- Unique terminal Slurm accounting: COMPLETED, exit 0:0, bizon, 180 CPUs,
  900 GiB; matching execution job, node, command, working directory,
  environment, source identity and finite ordered timestamps.
- The pinned dedicated Python runtime before/after admission and both
  checks recorded by the native runner; unchanged native Perl/system
  snapshots and prepared OrthoMCL sources.
- The retained BPO admission manifest, its completed 2-CPU/64-GiB job,
  frozen admission source, input records, index validation and unchanged
  failed-query diagnostics.
- Four unchanged staged BPO/index/GG copies in the native filenames, with
  original source hashes, single-link nonsymlink files and recorded mtimes;
  exact staging scope, validation command and output. These checks repeat
  after the independent index validation.
- One final-group output, unchanged complete native file inventory,
  nonempty graph/index/partition artifacts and the exact 3,003-file
  species-pair cache inventory. Inventory/hash checks repeat after the
  independent final-group and index validation.
- Complete final-group agreement with the native partition, GG ownership
  and input scope (984,137 proteins, 78 species), plus agreement with the
  runner's summary. Singleton omission is retained; graph edges are not
  substituted for final groups.
- An independent guarded native Storable index check on the staged copies.
  Native diagnostic output causes failure through the existing native-step
  helper. Source and input hashes are checked again before admission.

Checked records are deduplicated by path, rejecting conflicting identities,
to avoid repeatedly hashing identical large inputs within one verification
pass. This is not a reduction in the required artifact set.

Success is `corrected_orthomcl_native_outputs_admitted`, with pair semantics
`cross_species_final_group_cliques`. `accuracy_admitted` and
`publication_ready` remain false. Conversion and scoring are separate steps.
The report carries original query-failure evidence and does not imply that
conversion or clustering repaired missing BLAST hits. Shared-host native
time is not admitted as a matched end-to-end timing measurement.

## Verification

57 new tests cover execution and scheduler corruption, runtime identity,
staging provenance and mutation, hardlinks/symlinks, missing/extra/changed
native outputs, cache mismatch, conflicting records, pending-job rejection,
and refusal to overwrite an output path. Mocked complete orchestration
checks success, final-group mismatch and index mismatch, retaining failure
reports and source query diagnostics. These mocks do not claim a completed
production audit. Together with the 51 final-group tests, 108 focused tests
pass.

An actual dedicated-runtime invocation against pending native job 21750
passed runtime verification and rejected the missing COMPLETED accounting
before reading production outputs or creating the admission directory.
Batch script syntax passes `bash -n`.

Full unit suite with `ORTHOHMM_LEGACY_BLAST_SMOKE=1`: 3,967 passed in
86.37 seconds. Existing unrelated sample-output changes were not staged.

## Execution

Batch: `qfo_corrected_orthomcl_admit_batch_20260918.sh`, 2 CPUs/64 GiB/24 h,
no requeue, bizon. It runs a frozen validator under the dedicated Python
environment with isolated startup and a fresh bytecode-cache prefix.
The dependency must be `afterany:21750`, so a failed native job is observed
and rejected rather than silently accepted or resumed.

Expected output:
`benchmarks/work/qfo_corrected_orthomcl_admission_20260918/report.json`.
The submission identity and frozen validator revision are recorded below
after validation and submission. No corrected production admission or score
exists yet.
