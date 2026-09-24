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
