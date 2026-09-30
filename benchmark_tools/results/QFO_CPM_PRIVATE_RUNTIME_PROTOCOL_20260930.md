# High-CPM Private-Runtime Refinement Control

## Question

Does the unchanged original refinement-only CLI complete in the reconstructed
private Python/package environment independently validated by full OrthoBench
job 22377? The original admission 22155 and allocator diagnostic 22158 failed
during post-write collection/readback. Debugger 22159 and the forced-GC control
completed, but changed allocation or collection history. Do not repeat them.

This is a distinct runtime/environment arm, not a proposed scientific-method
change or another forced-GC control. It uses the original runner, including the
fresh `set(names)` after writing. No inserted stage logging, forced collection,
disabled GC, debugger, optimizer, new search or scoring is permitted.

## Fixed Inputs And Execution

Commit/push this protocol, controller and tests before one actual attempt.
`probe_cpm_private_runtime.py` requires this protocol's explicit SHA256 and
fixed hashes for the restored 22377 result, restoration plan, package audit,
private `pyvenv.cfg`, previous GC result and native refinement metadata.
Require the retained 22377 result to report independent reproduction equality.
Use its existing private inference interpreter, not a global package upgrade.

Validate the original runner SHA256
`a277ab01e7fbcbfaa15f7a63a092c78183baaabbbcc25cd23640f6a750eabbe2`,
frozen source modules, checkpoint/graph/seed inputs, earlier control records,
all retained reconstructed-base records and all matched installed package
payloads before/after execution. Canonicalize in-base interpreter aliases while
keeping their expected size/hash; reject escaped aliases and conflicting pins.
Generated wheel metadata, bytecode and relocated non-site wheel data remain
outside the matched package inventory, as in its retained validation.

Use a fresh output root with symlinks to unchanged seed/payload inputs. Run the
original CLI once with `--mode repeat-refinement` from the frozen replay launcher.
Inherit default GC behavior. A separate preliminary process must confirm Python
3.10.13, private NumPy 2.2.6, enabled GC and default thresholds 700/10/10; it does
not add imports or allocations to the refinement child's heap.

Child settings match allocator diagnostic 22158: `PYTHONMALLOC=debug`,
faulthandler, hash seed 0, no user site, frozen launcher on `PYTHONPATH`, and
one-thread OMP/OpenBLAS/MKL. Remove `PYTHONHOME`, `LD_PRELOAD`, `LD_LIBRARY_PATH`.
Use one permitted CPU, 64-GiB address-space ceiling, 300 CPU seconds, 360 wall
seconds (60 for the preliminary process), no core dumps. Timeout kills only the
owned child process group. The controller also confines itself to one CPU.
This direct shared-Threadripper diagnostic does not need the deferred timing
window; elapsed values are descriptive, never controlled resource comparisons.

## Prespecified Outcomes

Retain raw stdout/stderr, exit/signal/timeout, partial files and controller report.
There is at most one refinement attempt, no automatic retry or resume of a
failed attempt. Preflight failure records zero refinement attempts. Existing
output roots may not be overwritten.

Completion requires the original child JSON, exact native metadata apart from
the destination path, and independently read complete gene membership:
984,137 genes/memberships, 390,845 groups, 78 species, zero selected directed
refinement hits. Expected written partition is 23,875,927 bytes, SHA256
`f4c6f1973bc9636828baf1e6d9be3f416a18fc8ada5302fb502c082495fc1811`.
The controller's separate stdlib reader is additional output validation, not a
replacement for the original child's final reader.

Success is one observation in a different reconstructed environment. It cannot
prove memory safety, identify a CPython/package cause, repair original admission
22155 or release any downstream dependencies. Failure is also diagnostic only.
Both retain `seed_admitted=false`, `accuracy_evaluated=false` and
`publication_ready=false`. A separately reviewed admission contract remains
necessary before considering recovered high-CPM results for accuracy analyses.
No defaults, endpoints, scores, historical runtimes or completed jobs change.
