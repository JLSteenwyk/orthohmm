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

## Independent Native Output Admission

`admit_recovered_orthomcl.py` binds the completed native job and frozen runner,
revalidates recovered BPO evidence, and checks command/environment/timestamps,
runtime identities, staged copies, source reconstruction and final output/cache
inventories. It uses the existing independent final-group vs MCL partition/index
audit, then independently rechecks the native BPO indexes. Before success it
rechecks all provenance, staging, source/output inventories and runtimes.

Only `recovered_orthomcl_native_outputs_admitted` sets conversion authorization.
The pair semantics are explicitly cross-species final-group cliques; the report
does not claim these are native pairwise ortholog predictions. Search failures
remain intact and accuracy/publication flags remain false. Failed content or
post-validation checks leave a durable failed report, without conversion rights.

The isolated `qfo_blast_recovery_native_admit_20260923.sh` requires exact native
and BPO-admission executors plus the actual native job. It uses a separate
two-CPU/64-GiB allocation and fresh output/cache directories, and remains
unscheduled. The related contract/flow/content suites pass 205 tests, including
47 new tests. Shell syntax passes. Production validation and the subsequent
recovery-specific pair-conversion/scoring workflow remain outstanding.

## Recovered Pair Conversion

`prepare_recovered_orthomcl_pairs.py` now accepts only independently admitted
recovered native groups with conversion authorization. It binds the admission
job, source and executor, rechecks native scheduler accounting and all retained
records, and uses the pinned assessment reference mapping. A separate group
audit must reproduce the native admission's complete content summary before
the unchanged clique writer expands cross-species pairs. Pair totals and
mapping retention must agree exactly; no missing input hit is reconstructed.

The fresh output directory is `benchmarks/results/qfo_blast_recovery_pairs_v1`.
The successful state remains `recovered_orthomcl_pairs_prepared_unscored`.
Failed conversions retain reports and partial artifacts; no implicit retry or
replacement of the original pair workflow occurs. Query-failure coverage is
carried forward unchanged, and no accuracy/publication claim is authorized.

`qfo_blast_recovery_pairs_20260923.sh` supplies an isolated two-CPU/64-GiB,
no-requeue launcher with exact executor/job/report-digest arguments. It is not
scheduled. The focused suite passes 115 tests, including 25 new contract/flow
cases; shell syntax passes. Terminal conversion verification and a recovered
assessment/scoring handoff are still required before production QfO scoring.

## Frozen QfO Assessment Handoff

`run_qfo_recovered_orthomcl_assessment.py` now checks completed conversion
accounting and exact converter identity, pair counts/hashes and reference
retention, native admission and group-audit agreement, unchanged failed-query
coverage, and the pinned assessment environment/reference mapping. It calls
the existing assessment command builder without changing endpoints or settings.
Fresh paths are `qfo_blast_recovery_assessment_v1`, `qfo_benchmark/w/qc_mcr`,
and `qfo_benchmark/scoring/corrected_orthomcl_recovered`; participant identity
is `qfo_corrected_orthomcl_recovered`.

The eight-CPU/64-GiB no-requeue wrapper requires explicit conversion job,
manifest digest and exact clean executors. `--check-only` performs preflight
without starting scoring or creating output directories. A successful scoring
process is only `process_succeeded_pending_independent_admission`, with no
accuracy/publication authorization. The related suites pass 141 tests,
including 37 new contract/preflight/run cases; shell syntax passes.

No production scoring is scheduled. The recovery-specific independent score
admission adapter, actual completed upstream results, and final evidence
consolidation remain outstanding.

## Independent Recovered Score Admission

`admit_qfo_recovered_orthomcl_assessment.py` now validates the completed scoring
job, exact runner source and frozen executor, immutable preflight, rebuilt
conversion/command/reference provenance, complete output inventory and unique
native task trace. It uses the existing six-endpoint native metric validator
and rechecks input/output hashes before writing an admitted-score artifact.
No artifact is written on rejected validation; existing artifacts are not
overwritten. Admitted accuracy means validated benchmark scores, not independent
biological validation or publication readiness.

The assessment preflight can now reconstruct historical records from the frozen
scoring executor without demanding empty output directories. Launch still uses
the default fresh-directory requirement, and the historical options are not
exposed through its CLI. This logic is covered by an additional preflight test.
The adapter retains group/failure coverage, clique semantics and endpoint-specific
uncertainty caveats; the six-endpoint mean remains a secondary project summary.

`qfo_blast_recovery_score_admit_20260923.sh` requests two CPUs/64 GiB and requires
actual scoring/conversion jobs, exact executors and the pinned conversion digest.
The related suites report 97 passes, including 13 new adapter cases, with shell
syntax passing. It remains unscheduled. Production execution/admission, source
failure-impact analysis and final result consolidation are still outstanding.

## Frozen Executor

The common downstream code snapshot is
`benchmarks/work/publication_blast_recovery_downstream_v1_20260923` at
`a1dac0fd556833643476c77d5f9a5edaa0982fa2`.
[The receipt](qfo_blast_recovery_downstream_executor_20260923.json) records
eight source hashes, seven syntax-checked wrapper hashes and retained test/probe
artifacts. Sparse checkout status is clean. Use the original project root for
runtime/data paths, not the executor checkout as the data root.

Combined checkout tests passed 488 with two installed-fixture skips. The first
smoke-enabled attempt failed one checkout-local runtime lookup; its XML is
retained. A separate staged native probe executed from frozen code with the
real data root verified 12 groups and input preservation. Original data-root
installed staging/pair tests passed 53 with no skips. These fixture checks do
not admit production outputs. No downstream job is queued; parent jobs and
admission digests must be verified before each scheduled handoff.

## Frozen Consolidation And SwissTrees Executor

Use `benchmarks/work/publication_recovered_swiss_v1_20260923` at commit
`0b8efa9c6b20b5bdec41f798c136b54132e7db37` for comparison export and
`run_corrected_swiss_comparison.py` after independent recovered score admission.
This is separate from the native inference/scoring executor above. Its
[receipt](qfo_recovered_swiss_executor_20260923.json) records 123 passing tests,
source hashes and a real seven-method regression from the clean checkout.
All raw family counts, point estimates, comparisons and intervals exactly
match the retained seven-method analysis. No recovered scores were imputed.

Before a complete-panel run, verify the actual recovered score admission,
export a fresh comparison manifest using this executor, and pass that manifest's
actual SHA256 to the uncertainty runner. Keep the existing baseline and both
protocol files; the launcher checks their frozen hashes and all helper hashes.
Verify checkout commit and cleanliness again at execution. Use a fresh output
path, single-thread numerical-library settings and scheduled resources as in
the historical two-CPU/64-GiB SwissTrees batch. Do not run the historical batch
wrapper unchanged: it pins the old executor and seven-method manifest.

The primary stratified-analysis launcher remains tied to its historical sources
and three primary configurations; it is not an all-method recovery launcher.
Failure-impact analysis, production recovered scores, complete-panel uncertainty
and final publication consolidation are still outstanding.

## Failed-Query Membership Analysis

After whole-search and native-group admissions exist, use
`audit_orthomcl_failure_membership.py --search SEARCH_ADMISSION
--search-sha256 ACTUAL_SHA --native NATIVE_ADMISSION --native-sha256 ACTUAL_SHA
--output FRESH_OUTPUT` to join retained search diagnostics to the native
partition. It verifies matching query coverage, rehashes group-audit inputs,
reconstructs the full native partition and distinguishes final groups, MCL
singletons and proteins absent from the native index. It retains incoming,
outgoing and self-hit flags and does not invent singleton final groups.

This content audit relies on the supplied, pinned upstream admissions; it does
not rerun scheduler validation or the search. Its records can feed the separate
reference-exposure checker with `--allow-grouped-failed-queries` when grouped
failures exist. Mapping validity and reference exposure still require that
separate analysis. Neither incoming hits nor group membership establishes
repaired orthology or a counterfactual score. The module has fixture/CLI tests,
not a production recovered result.

Its frozen executor is now
`benchmarks/work/publication_failure_membership_v1_20260924` at
`92fd9e5ab6a5bbb2c048b2e4bbd87e25a9859253`. The
[receipt](qfo_failure_membership_executor_20260924.json) records source hashes
and 90 passing tests from that clean checkout. Recheck checkout identity and
cleanliness before use. Supply absolute paths to actual main-repository
admissions and a fresh output path; do not substitute partial batch reports.
Production execution remains pending complete search and native admissions.
