# Site and Frozen-Import Parser Controls

The default/debug no-site-import parser controls completed without reproducing
the crash. Next run exactly one `site_only` child and one `frozen_imports` child,
both with ordinary site initialization, debug allocator, hash seed zero,
faulthandler, no user site and one-thread numerical library settings. Use the
original replay cwd/PYTHONPATH and retained interpreter. No retry. Limit each
child to 120 wall/CPU seconds and 4 GiB address space; disable core dumps.

Both execute the same hash-pinned parser helper and saved partition. The first
adds site initialization but no explicit scientific imports. The second imports
the same replay, accuracy, refinement, partition audit and numeric checkpoint
audit modules as the failed refinement worker, requiring frozen source paths.
Do not load the numeric checkpoint, run refinement, construct a graph or call
an optimizer. Record resolved scientific source identities, all imported module
names, GC state/counters and before-import/before-parser/after-parser markers.

Check the previous control evidence, retained interpreter/libc, frozen source
files and recorded installed site-package identities before and after. Retain
each child exit, stdout/stderr, timeout or failure. A success only establishes
that this import/parser sequence completed once; it cannot exclude intermittent
memory corruption, reproduce the preceding full allocation history or admit
high-CPM scores. A failure is a localization observation, not automatic proof
that an imported module caused it. No scientific or frozen source is edited.
No DGX access or controlled timing is involved.
