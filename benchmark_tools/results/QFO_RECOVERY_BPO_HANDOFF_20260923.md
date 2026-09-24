# Recovery BPO Handoff

`prepare_blast_recovery_bpo.py` accepts only the recovery-specific admitted
search contract, not the interrupted original search or a candidate-only merge.
It requires the expected admission path and caller-pinned SHA256, completed
validator 22151, the frozen validator executor at
`a21c65f2b449d828e3fefc86148cbf4fd87b6ade`, unchanged transitive records,
the admitted merge identity and exact corrected database parity.

The wrapper uses the existing dedicated Python runtime verifier and unchanged
`prepare_orthomcl_bpo_checkpoint.py`. It does not modify the converter,
1e-5 E-value cutoff, native Perl indexing or independent index/content checks.
Inputs remain the recovered candidate table and original corrected all.fa.
Output is a fresh `benchmarks/results/qfo_blast_recovery_bpo_v1` directory;
the old partial search and old BPO workflow are not overwritten or relabeled.

Source HSP and directed-pair totals must match the admitted search. Runtime
and input identities are rechecked after checkpoint preparation. Successful
preparation still requires independent terminal checkpoint admission before
native inference; accuracy, publication and downstream authorization remain
false. Logged query failures remain in the provenance and are not repaired
by conversion.

All 122 focused recovery tests pass, including fifteen new input-contract
and orchestration cases. These fixture tests cover both successful preparation
and checkpoint/count/runtime failures, but do not replace production native
conversion or validate a real future report. The new wrapper is not scheduled.
It requires a real successful search-admission report and pinned digest, plus
frozen scheduling and an independent recovery checkpoint admission workflow.
Existing held downstream jobs remain untouched.

## Isolated Launcher Verification

The recovery wrapper now follows the original script's explicit repository
import-path setup, allowing direct execution with `python -I -B` from outside
the repository. Its isolated `--help` entry point was tested. A production-style
runtime probe using the dedicated Python 3.10.13 environment, minimal environment
variables and an unused `-X pycache_prefix` verified all 2,928 runtime records.
An initial probe without that cache prefix was correctly rejected; the runtime
contract was not relaxed.

`qfo_blast_recovery_bpo_prepare_20260923.sh` mirrors those verified launcher
settings and requires an explicitly supplied future admission digest. Shell
syntax passes. Additional scheduler/executor/source rejection tests bring the
focused recovery suite to 128 passing tests. This launcher remains unscheduled;
the runtime probe did not convert any production BLAST rows or admit a checkpoint.

## Independent Recovered Checkpoint Admission

`admit_blast_recovery_bpo.py` now implements the recovery-specific independent
checkpoint gate. Its caller must supply the actual completed preparation job,
clean retained executor and exact revision. The gate pins preparer and reused
validator source hashes, binds the recovered search admission and query-failure
coverage, verifies checkpoint inventory and native commands, and rechecks the
full BPO content plus native Perl indexes in a fresh directory. The dedicated
Python, Perl and system-helper identities are checked before/after validation.

Only successful validation yields `checkpoint_admitted=true`; accuracy,
publication readiness and automatic downstream authorization remain false.
Failures after validation starts produce a durable failed report. No original
held job is released. The entry point supports isolated direct execution.

The focused preparation/admission/content/index suite has **123 passing tests
and four skips** for unavailable opt-in native integration dependencies. The
new contract tests cover failed/pending preparation, wrong allocation, source,
checkpoint, runtime, timestamps, and premature accuracy/downstream flags.
These checks do not replace real-data execution or full recovered-admission
orchestration testing. No independent admission wrapper has been scheduled;
the future preparation job, executor and result do not yet exist. Next work is
orchestration failure-path coverage, then freezing the actual preparation and
admission executors after search admission. Native inference still needs a
recovery-aware handoff and separate validation.

## Flow And Native Fixture Checks

Added twelve recovered-admission orchestration cases using real temporary file
hashes and mocked expensive scientific/runtime operations. They verify success,
refusal to overwrite, pending/wrong-revision gates, checkpoint/search/coverage
and source-count failures, missing search provenance, failed content recheck,
changed input bytes, and post-validation Python/native runtime failure. No
failed case admits the checkpoint or authorizes downstream execution; successful
fixtures retain failed-query accounting.

The new flow/contract/preparation subset passes **47 tests**. Separately enabled
`ORTHOHMM_LEGACY_BLAST_SMOKE=1` for the existing native index, preparation and
independent-check suites: **68 passed with no skips**, including installed
OrthoMCL/BioPerl fixture indexing and revalidation. This resolves the earlier
opt-in test gap; it does not validate a future production checkpoint.

The new `qfo_blast_recovery_bpo_admit_20260923.sh` requests two CPUs, 64 GiB,
24 hours and no requeue, using the pinned dedicated Python with `env -i`,
isolated imports and a fresh bytecode-cache prefix. It requires explicit future
preparation job/executor/revision, checks its own frozen executor, and writes
only to a fresh recovery-specific admission directory. Shell syntax passes.
It remains unscheduled; the required production inputs are not yet available.

## Native Input Evidence Gate

`recovered_orthomcl_inputs.py` verifies the recovery-specific checkpoint report,
explicit successful admission/preparation jobs, pinned admission source and
clean executor, full recovered search provenance and unchanged failure counts.
BPO/index paths must match the recovered checkpoint, with consistent nonempty
records. The species mapping comes from the pinned original prepared-input
manifest: it is not assumed to be present in BLAST recovery's record list.
The actual mapping checksum and 984137-protein/78-species scope were verified.

The gate's 27 tests and the existing preparation/admission tests total 74 passes.
This is an evidence-only component, not a native launcher. It authorizes no
execution and does not overwrite original outputs or release old jobs. Before
launch, connect dedicated runtime checks and fresh staging to the recovered
evidence, and add an explicit admission-job identity to the admission report
(currently the job is supplied by the caller and verified separately in Slurm).
Native final-group admission remains a separate required step.

## Admission Job Identity Bound

The admission report now embeds its own Slurm job ID, node, CPU allocation and
memory allocation, distinct from the completed preparation job's accounting.
Admission refuses unscheduled or incorrectly allocated execution. The isolated
wrapper explicitly preserves those Slurm variables through `env -i`.
The native input gate requires the report's identity to match the explicit
completed admission job and pins the revised admission source hash. This closes
the job-identity limitation recorded in the preceding section; no production
admission report existed before this change.

Eight execution-identity cases and four downstream mismatch cases extend the
focused suite to **86 passing tests**. Both successful and failed validation
reports retain their own job ID. Shell syntax also passes. The launcher/runtime
integration and final-group validation for native inference remain outstanding;
this change does not schedule admission or authorize native inference.

## Separate Native Launcher

`run_recovered_orthomcl.py` now connects the evidence gate to dedicated runtime
verification, the existing source configurator, fresh input staging, native
mode-4 inference, and existing final-group/cache checks. It writes only under
`benchmarks/results/qfo_blast_recovery_native_v1`, refuses existing directories
or symlinks, and retains failure reports without implicit retries. Native
sources come from the checksum-pinned original manifest; only fresh local
paths are configured. Scientific defaults and the 64-worker pair patch are
unchanged. Input hashes and mtimes, runtime identities and expected caches
are checked after inference.

`qfo_blast_recovery_native_20260923.sh` requests 180 CPUs/900 GiB and no requeue.
It requires the actual completed BPO admission job, explicitly pinned report
digest, and exact clean native/admission executors. It remains unscheduled.
Success is only `recovered_native_exited_zero_pending_admission`; independent
terminal/final-group validation and scoring remain necessary. The old held
inference job is not reused or released.

Recovery-focused tests: 105 passed, including 19 new orchestration/preflight
cases. Reused native staging/configuration/original-runner/pair-parallel tests:
71 passed with native smoke enabled, no skips. Shell syntax passes. These
fixture results do not establish production inference success or accuracy.
