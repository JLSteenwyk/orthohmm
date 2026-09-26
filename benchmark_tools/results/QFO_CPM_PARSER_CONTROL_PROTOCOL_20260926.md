# Parser-Only High-CPM Controls

Before executing either arm: test the diagnostic and freeze this protocol in
its checked input records. Run exactly one fresh process with the default
allocator and one with `PYTHONMALLOC=debug`, both with `-B -S`, hash seed zero
and faulthandler enabled. Each has a 120-second wall/CPU limit, a 4-GiB address
space limit and no core dump. No retries. No optimizer or refinement executes.

Use the retained interpreter/libc identities, exact frozen `read_partition`
function, gene universe and refined output partition. Extract only the function
AST from the hash-pinned local module; do not execute its imports or other
module-level code. Leave garbage collection enabled at its existing thresholds.
Read the partition once, retain GC counters, module inventory, coverage and
all stdout/stderr, including nonzero exits/timeouts. Preserve failed admissions.

This targets the Python parsing/allocation boundary identified by the earlier
faulthandler stack. It deliberately omits earlier allocation history, scientific
imports, array conversion and native refinement. A successful read cannot
exclude native corruption or identify the earlier crash cause; a failing read
would motivate a parser/interpreter-focused reproducer. Neither result admits
high-CPM scores, changes defaults or releases dependent jobs. Shared-host
elapsed times are descriptive. The DGX remains deferred.
