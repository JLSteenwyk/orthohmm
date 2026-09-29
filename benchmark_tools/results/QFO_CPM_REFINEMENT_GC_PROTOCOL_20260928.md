# Refinement Garbage-Collection Boundary Diagnostic

## Question

Does the retained high-CPM refinement state fail during a forced generation-2
collection before refinement, after refinement returns, or after its output
is written? Existing isolated parsing, scientific-import and checkpoint/setup
controls completed; the failed admission and allocator diagnostic crashed
during later parser-triggered collection. This tests the unisolated refinement
boundary, not a proposed fix or a new accuracy evaluation.

## Fixed Execution

Run `probe_cpm_refinement_gc.py` once after committing and pushing this protocol,
runner and tests. Use the retained interpreter, frozen replay imports, numeric
checkpoint, high-CPM seed partition and graph arrays. Validate the earlier
checkpoint/setup receipt and all its checked records before/after execution.
Require the protocol digest explicitly. Do not execute Leiden or rerun search.

Child settings: debug allocator and faulthandler as in the prior setup control,
one CPU selected from the inherited permitted affinity, one-thread numerical
libraries, 64-GiB address-space ceiling, 300 CPU seconds, 360 wall seconds,
disabled core dumps. GC remains enabled at its original 700/10/10 thresholds.
This is a bounded direct subprocess diagnostic on the shared Threadripper,
not a Slurm production timing identity or a controlled performance measurement.

Load and validate the checkpoint and seed groups; map graph arrays; select
the frozen refinement hits. Force collection immediately before refinement,
immediately after it returns and immediately after writing the new partition.
Record flushed before/after markers and GC statistics around every collection,
refinement, write and final readback. Keep the actual scientific objects alive
through readback. Do not disable GC, use a different reader, or change scores,
thresholds or refinement settings.

## Outcomes

Retain stdout, stderr, return code, timestamps and any partial output regardless
of success, signal or timeout. There is one child attempt and no retry. A
successful result must cover 984,137 genes in 390,845 groups, 78 species and
zero selected directed refinement hits; its written partition must match the
retained refinement SHA256
`f4c6f1973bc9636828baf1e6d9be3f416a18fc8ada5302fb502c082495fc1811`.

A failure localizes the last observed stage only; it does not identify the
corrupting operation. A success does not exclude intermittent corruption or
reproduce the original allocation history: forced collections and logging
deliberately perturb that history. Neither outcome repairs or supersedes failed
admission 22155, admits high-CPM scores, or authorizes a retry. A native crash
trace or a separately justified narrower diagnostic would still be needed.
