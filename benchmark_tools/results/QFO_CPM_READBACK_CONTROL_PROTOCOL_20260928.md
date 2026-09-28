# Isolated Partition Readback Control

One prospective diagnostic, not a scientific retry or admission change.
The retained GC-stack SIGSEGV occurred during post-write `read_partition`.
GDB did not reproduce it, and startup Memcheck reports also occurred without
scientific imports. Test whether the saved reader and saved partition bytes
alone reproduce the failure, without refinement or scientific extensions.

Extract only the `read_partition` function AST from its checksum-pinned frozen
module; retain the function body and source line numbers, omit module-level
imports. Use the retained Python executable with `-B -S`, allocator debug,
faulthandler, enabled/default garbage collection, and the original replay cwd
and PYTHONPATH. Verify gene names, initial partition, failed run's refined
partition, interpreter and source identities before/after. Read the initial
partition, discard it, then read the refined partition and retain the groups.
Require 984,137 unique genes and 390,845 refined groups. Record GC statistics
and loaded modules. Do not force or disable GC, optimize, refine, or score.

Execute once, no retries, under local Slurm: one CPU, 8 GiB, ten-minute
allocation, five-minute child timeout, no requeue. Freeze and push protocol,
runner and focused tests before submission. Preserve nonzero exits/timeouts
and incomplete evidence. Do not overwrite prior artifacts.

A crash is evidence for a smaller reproducer, not proof of its cause. A
successful read means only that this isolated allocation/import history did
not reproduce the crash. It does not test the preceding refinement heap state,
exact GC schedule, native extensions or complete runtime. File-byte identity
does not establish runtime-library equivalence. No success rate, safety,
repair, accuracy, controlled timing or downstream admission claim follows.
Retain failed admission 22155 and missing high-CPM endpoint regardless.
