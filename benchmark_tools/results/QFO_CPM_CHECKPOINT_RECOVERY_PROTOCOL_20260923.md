# High-CPM Checkpoint Recovery Protocol

## Scope and Decision

Prespecified before any recovered high-CPM partition or accuracy score is
observed. Permit one checkpoint continuation, conditional on the admission
gates below, of failed job **22081_1**. This is not permission to relabel that
job successful, repeat trials until a favorable score appears, or change the
scientific method. Preserve every existing artifact and the failed scheduler
record. DGX work remains deferred.

The earlier boundary protocol stopped before optimization. Diagnostic **22152**
subsequently reached that boundary with exact saved graph parity; it did not
run Leiden or authorize recovery itself. This separate protocol defines the
proposed scientific continuation and its required validation.

## Source Review

Reviewed the frozen launcher
`benchmarks/work/publication_qfo_replay_native_v1/benchmark_tools/replay_high_sensitivity.py`,
its CPM context, the completed-stage auditor, and candidate preparation code.
The high-CPM output root is
`benchmarks/results/qfo_parameter_cpm_replay_v3/cpm_high`.

The saved worker reports stages 0 initial, 1 multipass and 2 profile_base as
checked; stage 3 profile_expanded failed with SIGSEGV. Its complete input graph
was saved before worker execution. The HMM expansion and final singleton-edge
construction therefore precede the failed call. The two multipass output files
exist, each 23,875,927 bytes. No final `replay.json` exists, and neither profile
partition has been written. These are filesystem/control-flow observations,
not fresh admission of any retained stage.

After stage 3, the frozen source reads the partition, writes strict_profiles,
then calls `refine_cluster_indices` using the original numeric hit checkpoint,
production refinement hits and the final graph arrays. With one profile
iteration, there is no subsequent HMM search. Thus these saved inputs can
support a scientifically equivalent continuation without reconstructing the
in-memory profile objects. They do not recover lost profile counters or timings.

The existing admission requires successful job 22081_1 and a complete original
replay report. The candidate driver also pins that admission and its seed path.
Neither contract can accept a recovery by changing only a job dependency.

## Required Preflight

1. Recheck the corrected replay plan, numeric checkpoint, gene order, all input
   FASTAs, frozen scientific imports and native runtime using existing gates.
   Require the same 984,137 genes, 78 species, CPM 0.12, seed 4, isolates enabled,
   BLOSUM62, one profile iteration and profile_min_species 1.
2. Independently revalidate stages 0 through 2, including original/copied payload
   hashes, metadata, commands, execution records, native graph observations,
   checked partitions and full unique gene coverage. Compare execution records
   with the failed parent worker. Do not require or invent successful stage 3.
3. Recheck the failed stage's saved payload and original manifest, the failure
   report, and boundary result 22152. Require 25,501,180 edges and the already
   recorded complete endpoint and weight fingerprints. Verify all source,
   adapter and native library identities, not just graph dimensions.
4. Reconstruct production refinement hits from the admitted numeric checkpoint.
   Recompute the existing multipass-refined partition with the frozen function
   and compare memberships exactly to the retained output. This tests the
   continuation's refinement input and call conventions before recovery.
5. Use a new output root and a clean, committed detached executor. Bind the
   protocol and fresh preflight report by digest. Reject existing output paths,
   changed inputs, missing evidence or any failed gate. Do not overwrite the
   original graph, metadata, partitions, reports or executor.

## Execution and Validation

Use the frozen Leiden worker with the previously checked Python-pair constructor
adapter, one-CPU affinity and thread limits of one. Only redirect output metadata;
leave graph order, weights, seed and resolution unchanged. Allocate one CPU,
64 GiB, 24 hours, no requeue, on the existing host. Record this as incremental
shared-host recovery, never end-to-end controlled timing. If measured preflight
memory requires a larger allocation, document an amendment before submission.

Allow the single stage-3 optimization to complete. Capture the actual graph
handed to Leiden and require exhaustive endpoint/weight parity with the saved
input, alongside native runtime and partition checks. A preoptimizer diagnostic
alone is not sufficient. Retain failure and stop on any native error or mismatch;
do not automatically resubmit, vary parameters or select among outcomes.

Read the recovered partition, write strict_profiles, then execute the exact
frozen refinement call from the source review. Require full, duplicate-free gene
coverage for both partitions. Independently repeat refinement from the recovered
partition and numeric inputs and compare memberships. Retain both original
multipass outputs with their actual provenance. Report stages as reused or newly
computed rather than pretending all four ran in one successful invocation.

The recovery report must have a distinct schema/status and identify the original
failure, preflight, scheduler job, executor, commands, settings, all input/output
hashes and measured incremental resources. Unavailable HMM profile counters,
original successful-stage timing totals and full-run totals remain explicitly
unavailable; do not estimate them or synthesize an original `replay.json`.

## Admission and Downstream Work

A separate process must revalidate the completed recovery scheduler record,
preflight evidence, original stages, actual optimizer input, recovered partition,
refinement reconstruction and before/after hashes. Recovery completion alone
does not admit a candidate seed or accuracy result.

Implement explicit recovery-aware candidate admission and subsequent provenance
handoffs. Reuse the frozen candidate expansion, phylogeny, pair conversion and
scoring logic without changing scientific parameters. Do not relax original
validators, impersonate their statuses or release the blocked high-arm chain
against the old failed job. Retain low-CPM and the original failure untouched.

Only after independent recovery and downstream admissions may high-CPM enter
the existing seven-arm parameter table and paired family uncertainty analysis.
Keep the prespecified multiplicity family, missing-result handling, endpoints
and seeds unchanged. Report the failure and recovery alongside any final score;
no accuracy result from this development-exposed arm establishes generalization.

## Implementation Gates

Before submission, test altered payloads/settings/runtime, missing or changed
predecessor partitions, native failure, absent graph observations, partial or
duplicate-gene outputs, refinement mismatch, accidental overwrite and scheduler
failure. Include a small real frozen-worker optimizer/refinement integration
test, not only mocked observations. Commit and push the tested implementation
and this protocol before freezing the executor. No recovery job was submitted
when this protocol was written.
