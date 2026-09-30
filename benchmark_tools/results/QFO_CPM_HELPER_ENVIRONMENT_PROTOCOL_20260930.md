# High-CPM Helper-Complete Environment Arm

## Rationale And Frozen Selection

The preceding private-runtime CLI failed at the frozen benchmark helper's
Biopython import, before refinement. Preserve that failed attempt and the
[independent readback](qfo_cpm_private_runtime_readback_20260930.json).
Do not rerun it unchanged or add packages to the validated full-OrthoBench
inference environment. This new arm prepares a separate helper-complete venv
before one unchanged-runner attempt; it is not a scientific-method change.

Use the reconstructed CPython 3.10.13 base and exactly the eleven inference
wheels already bound by restored OrthoBench 22377's plan, plus its retained
reader Biopython 1.87 wheel. All twelve are local, already acquired and
hash-pinned. No new network acquisition, source build, global installation,
version optimization or dependency resolution. The original inference wheels
retain NumPy 2.2.6, Numba 0.67.0, llvmlite 0.49.0, DendroPy 5.1.0,
igraph 1.0.0, Leiden 0.11.0, python-igraph 1.0.0, texttable 1.7.0,
OrthoHMM 0.5.0, pip 26.2.1 and setuptools 83.0.0.

Biopython's wheel SHA256 is
`1e951f4862ffc1dccc28e1c25245059fa653d86028a88f5dfe1b7875875f3a4f`.
Selection checks both retained wheel metadata and expected byte identities.
No inference/reader package may be silently substituted or combined through
`PYTHONPATH` or system site packages. Scientific imports must remain in the
frozen QfO replay source, not the installed OrthoHMM wheel.

## Preparation Gates

Commit/push `prepare_cpm_helper_environment.py`, amended controller, tests and
this protocol before actual preparation. Require this protocol's explicit SHA256
and the fixed prerequisite hashes already in the controller. Use a fresh output
root, copied interpreter and `include-system-site-packages = false`.
Install staged local wheels offline with hashes, binary-only, no dependencies
and no bytecode compilation; run `pip check`.

In a separate preliminary process, import every module used by the original
worker's import block, including replay, historical reader, stage diagnostics,
accuracy, refinement and NumPy. Those imports pull their actual transitive
imports. Require the exact twelve-package inventory, reconstructed base prefix,
enabled GC with thresholds 700/10/10 and no observed module outside the new venv,
reconstructed base or frozen replay source. Bind every observed module file.
This observes this import path, not every optional branch or dynamic-library
dependency. No checkpoint loading, refinement, optimization or scoring yet.

Compare installed package payloads to all twelve wheels. Preserve generated
RECORD/bytecode and relocated non-site payload exclusions. Check the complete
retained base inventory, source/helper/prerequisite records, copied wheel hashes,
interpreter identity, venv configuration, import report and installed payloads.
Retain preparation logs/failures; existing roots cannot be overwritten.
Children use one CPU, 64-GiB address space, 300 CPU seconds and 180 wall seconds
per preparation stage, debug allocator/faulthandler, one-thread numeric libraries,
no user-site/global environment injection. No automatic preparation retry.

## One Refinement Attempt

After successful preparation and independent stdlib readback, commit/push the
preparation result before execution. Invoke `probe_cpm_private_runtime.py` with
this protocol digest plus the new preparation receipt path and exact digest.
The controller rechecks that receipt, its inputs/payloads and import report,
and binds the new interpreter to the same reconstructed base. It also retains
the existing original source/data/prerequisite checks. There is at most one
native refinement CLI attempt in a new output root, with no automatic retry.

All original scientific arguments, frozen runner/source/checkpoint/graph/seed
bytes and its fresh post-write `set(names)` remain unchanged. Do not instrument
the worker, force/disable GC, substitute a reader, run a debugger or optimize
again. The preliminary import probe is not executed inside that child's heap.
Use the controller's existing one-CPU/64-GiB/300-CPU-second/360-wall-second bounds
and exact metadata/partition/full-membership completion checks from the
[original protocol](QFO_CPM_PRIVATE_RUNTIME_PROTOCOL_20260930.md).

Retain all logs, return/signal/timeout, partial files and result regardless of
outcome. No outcome identifies the cause of the historical crash or repairs
failed admission 22155. All seed/accuracy/publication admission flags stay false;
a separate prospective admission contract is still required for recovered
accuracy results. No historical result, score, endpoint, default, timing identity
or dependent job changes. Shared-host elapsed values are descriptive only;
controlled timing remains deferred and no DGX or unrelated-process action occurs.
