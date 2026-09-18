# Corrected HMM Native Admission Queued

Job **21720** is pending `afterany:21706_0` on bizon, with 2 CPUs,
64 GiB and a four-hour limit, without requeue. It validates inference's
numeric scheduler job **21707**. The array dependency was confirmed with
`scontrol show job 21720`; no inference was restarted or changed.

The detached validator executor is
`benchmarks/work/publication_qfo_corrected_hmm_admission_v1`, commit
`7b5214a5d2169338c4cbe80293fae491ee5b951f`.
Batch: `qfo_corrected_high_sensitivity_admit_batch_20260918.sh`.
Expected output:
`benchmarks/work/qfo_corrected_high_sensitivity_admission_20260918.json`.
Log: `benchmarks/work/qfo_corrected_hmm_admit_21720.log`.
The output does not yet exist; nothing is admitted by this submission.

## Validation Scope

The validator requires successful terminal scheduler status and the exact
32-CPU bizon inference identity, frozen primary plan, original inference
executor, command/cwd, timestamps and full native output inventory. It
reruns the existing input/runtime verifier, then checks complete corrected
gene/species coverage, checkpoint integrity and raw-cluster equivalence
of exported groups (allowing the native singleton completion).

Native metrics live beside, not inside, the output directory. Their content
is hashed separately and their command/cwd/timestamps must agree with the
execution. This post-completion hash is not a hash recorded by the native
runner. The module invocation is checked against its `__main__.py` metrics
representation without discarding the remaining command arguments.

During review, the initial content validator's synthetic fixture was found
to flatten native settings incorrectly. Both frozen and current metrics
writers put settings under `metadata`. The validator and fixture now use
the actual schema, with an integration test invoking `PipelineMetrics` and
a negative test rejecting flattened settings. No corrected output had been
admitted using the initial validator.

All 76 focused admission/content/checkpoint/ownership tests pass. The real
historical QfO metrics pass the corrected settings/count check (976,504 genes,
78 species, 390,817 groups). Direct invocation against still-running 21707
rejected it before reading partial artifacts and created no report. These
checks do not substitute for running the entire admission on completed
corrected output, or establish biological accuracy/search completeness.

After successful admission, freeze the corrected factorial/satellite replay
against the actual checkpoint manifest hash. Conversion and scoring need
their own validation. Shared-host measurements remain descriptive; dedicated
DGX timing admission is separate.
