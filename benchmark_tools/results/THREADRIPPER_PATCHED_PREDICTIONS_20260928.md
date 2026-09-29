# Patched Native Fixture Predictions

The [successful receipt](threadripper_patched_predictions_20260928.json) records
fresh native inference for both OrthoHMM timing configurations using the separate
patched interpreter. It preserves the first two frozen command-plan settings,
including 32 requested CPUs and four threads per worker, changing interpreter,
fixture input and output locations only. No full dataset or production timing
identity was run. The controller ran outside Slurm on the shared host; these
executions provide no controlled timing evidence.

| Configuration | Genes | Groups | Native phylogenetic pairs | Byte-identical prediction files |
|---|---:|---:|---:|---:|
| High sensitivity | 16 | 3 | Not applicable | 1 |
| satellite_v2 with inferred phylogeny | 16 | 3 | 36 | 3 |

Comparisons use retained outputs from fixtures 22367/22368, not freshly selected
reference answers. High-sensitivity orthogroups match; phylogenetic orthogroups,
root groups and pair-table bytes match. Validation additionally checks exact
input coverage, native command/cwd, completed metrics, group counts and sorted
unique native pairs. Forty frozen source records, prior evidence, candidate
receipt, plan and interpreter identity are checked around the executions.
Native MAFFT/FastTree/DIAMOND entrypoints are checked before invocation.

## Retained Preparation Failures

The [first attempt](threadripper_patched_predictions_initial_failure_20260928.json)
failed before native launch because the driver passed a persistent input path
to the existing helper that requires tmpfs. Its
[executed driver](probe_patched_timing_predictions_initial_20260928.py) is retained.

The [second attempt](threadripper_patched_predictions_v2_failure_20260928.json)
launched the CLI, which printed `Output directory does not exist` and exited
zero before inference. The output validator rejected the missing metrics.
Its [executed driver](probe_patched_timing_predictions_v2_20260928.py), local
native log and command receipts remain intact. The corrected driver creates
the required fresh output directory before launch. The third attempt then ran
both methods through inference once; no completed inference result was retried
or discarded to obtain a match. Frozen scientific code was not changed to
alter the CLI's historical exit behavior.

## Limits And Next Action

Nineteen focused tests pass across prediction comparison, candidate building
and import comparison. Tests include changed prior evidence, differing output
bytes, empty comparisons and phylogeny subdirectory handling; they do not
constitute exhaustive orchestration coverage.

This small fixture does not exercise satellite merging, full-scale memory
behavior or every optional dependency branch. It does not establish complete
runtime isolation, continuous process-tree containment, current whole-environment
security or genome-scale equivalence. No OrthoFinder rerun was needed for this
OrthoHMM interpreter-only parity check. Next, bind the new deployment to the
collector and validate collector v5 on native fixtures; final resource accounting
and verified quiet-host eligibility remain required before the 27 timing runs.
