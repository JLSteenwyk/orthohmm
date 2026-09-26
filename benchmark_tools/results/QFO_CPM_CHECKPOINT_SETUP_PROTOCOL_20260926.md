# One Checkpoint and Pre-Refinement Setup Control

After successful parser-only and import-context controls, run exactly one
debug-allocator child with the same frozen replay imports, cwd, interpreter,
hash seed and one-thread settings. Bound it to 300 wall/CPU seconds, 8 GiB
address space and no core dump. No retries and no optimizer or refinement.

Recheck preceding control records and the frozen source/runtime identities,
including all five scientific modules from the last child. Bind the numeric
checkpoint to its pinned plan and manifest and check all five checkpoint files.
Bind the saved seed partition and graph arrays to exact checksums. Preserve
the source payloads read-only.

Within the child, call the original audited numeric checkpoint loader. Record
a flushed boundary marker, then parse the retained refined partition using the
original frozen parser. Next validate gene order, parse/index the seed groups,
memory-map the graph arrays and call `production_refinement_hits`. Keep these
objects alive while parsing the retained refined partition again. Stop before
`refine_cluster_indices`. Report shapes/dtypes/mapping status, coverage, hit
counts, original checkpoint evidence, scientific source identities and GC
state/counts at each marker. Do not change GC settings or force collections.

The extra parser probes and diagnostic state alter allocation history; even
success does not reproduce the failed refinement process, establish memory
safety or identify its cause. Failure should retain the last boundary and
signal/timeout without automatic resubmission. Neither outcome repairs failed
admission 22155 or releases candidate 22156. Shared-host times are descriptive.
No default changes, accuracy scores or DGX access are involved.
