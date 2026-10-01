# Recovered High-CPM Phylogeny: Runtime Gate Rejected

The separate recovered-candidate handoff is implemented and tested, but **no
phylogeny job was submitted**. A live read-only `verify_sources(Path.cwd())`
preflight on October 1 rejected the full scientific runtime before creating
the new output directory or invoking fresh candidate admission/native inference.
This is an unscheduled preparation failure, not a scientific benchmark attempt.

The [machine-readable receipt](qfo_cpm_helper_phylogeny_preflight_20261001.json)
binds the actual command, scheduler observation, input/source/protocol identities,
metadata differences and absent output. Final focused validation: 199 tests
passed in 3.39 seconds, including the actual fail-closed inventory regression;
shell syntax passed. This is not a rerun of the full repository test suite.

## Actual Check

The command used `/home/bizon/anaconda3/bin/python -B` from the original
repository root, with PYTHONNOUSERSITE/PYTHONPATH/PYTHONHOME/LD_PRELOAD/
LD_LIBRARY_PATH unset and hash seed 0, OMP/OPENBLAS/MKL thread limits 1.
It invoked `benchmark_tools.run_helper_cpm_phylogeny.verify_sources` directly;
no Slurm task, output argument, optimizer or scientific stage was executed.

Current sacct confirms candidate construction 22385 and admission 22386 both
COMPLETED 0:0 with 2 CPUs/64 GiB on bizon. Original 22155 remains FAILED 1:0;
22156 remains CANCELLED with zero execution. The recovered candidate readback,
full report, source executor/Git bindings and bound records passed before
`run_qfo_parameter_phylogeny.verify_baseline` reached `verify_environment`.
That unchanged strict check raised:

```text
ValueError: Package inventory changed: orthohmm
```

A separate read-only importlib.metadata inventory comparison found identical
Python version text, no removed distributions, 607 expected versus 611 observed
distributions, four additions and one changed version:

| Distribution | Frozen Metadata | Current Metadata |
| --- | --- | --- |
| ProDy | absent | 2.6.1 |
| dm-tree | absent | 0.1.10 |
| ml_collections | absent | 1.1.0 |
| propka | absent | 3.5.1 |
| pyparsing | 3.2.1 | 3.1.1 |

All five observed metadata paths are in the shared
`/home/bizon/anaconda3/lib/python3.10/site-packages`, not just optional
user-site packages. No shared installation was modified. No inference accuracy
or runtime effect follows from this inventory difference alone.

This is the same already documented
[shared dependency drift](THREADRIPPER_DEPENDENCY_DRIFT_20260928.md).
That earlier audit, not a new exhaustive scan here, found changed historical
pyparsing import bytes. Reuse its retained receipt; do not repeat the broad
runtime scan. The new preflight is evidence the original full-runtime gate
remains unsatisfied, not a reason to weaken it or hide the changed inventory.

## Next Action

Retain the tested strict entry point and prospective protocol. No launch under
that protocol is admissible in the current shared environment. Integrate a
separate explicit **full inferred-phylogeny deployment amendment** using the
already validated private runtime work, with complete source/runtime/input
bindings and actual native output parity evidence. The helper-complete runtime
amendment remains limited to independent refinement and does not authorize
downstream inference. Do not silently substitute it.

Reuse [private runtime evidence](THREADRIPPER_PRIVATE_RUNTIME_20260928.md),
its linked later packaging/deployment updates and the
[completed archive reproduction](RESTORED_ARCHIVE_FULL_OB_RESULT_22377.md)
where bindings and scope apply. Version labels or an unrelated reproduction
alone do not establish the required full QfO handoff. Preserve the original
full-baseline verifier, failed/cancelled gates and all scientific settings.

Next require a hash-bound private-environment handoff, launcher/import identity,
raw-checkpoint compatibility and fixed-configuration native parity before one
new local phylogeny attempt. Candidate content remains admitted but unscored;
native output admission, pair conversion/scoring and the fixed seven-arm/
18-endpoint analysis remain open. Controlled timing is deferred. No quiet-host
question, contention probe, DGX access or unrelated job/service action occurred.
