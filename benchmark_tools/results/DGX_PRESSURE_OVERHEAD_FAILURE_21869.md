# Pressure Panel Deployment Failure

All18tasks of21869 terminated FAILED1:0 after4-5scheduler seconds.
Accounting contains no native srun measurement steps. Recorder21870
completed0:0, retaining all18detailed terminal records across39polls with
zero observation errors. The entire panel was allowed to terminate in its
frozen order; no selective cancellation, restart or DGX inspection occurred
during the quiet window. That window is now closed.

## Diagnosis

After terminal accounting, collected all batch logs and run directories.
Each `verification.json` reports `verified_wrapper_failed`,
`ModuleNotFoundError`, `No module named 'numpy'`. No run has a preparation
record, measurement directory or native measurement result. Thus0/18native
measurements exist; no overhead ratio or equivalent-work comparison can be
computed from this panel.

The submitted batch script used `/usr/bin/python3`. Its earlier non-executing
preflight selected tasks and verified recipe bytes but did not invoke the
frozen input enumerator. That enumerator imports NumPy, which is absent from
the DGX system interpreter. The script now uses the already pinned
`/home/jlsteenwyk/projects/orthohmm-publication/envs/orthohmm/bin/python`.
No package was installed and no runtime manifest or environment was changed.

A new [non-executing enumeration check](dgx_pressure_enumerator_preflight_20260919.json)
using that exact interpreter succeeds and reproduces the original four-file
input order. It reports Python3.10.13. The batch logs also contain Slurm's
fallback from an unavailable inherited TMPDIR to `/tmp`; this warning is
retained and is not the recorded Python failure cause.

## Evidence and Repair

[Failure summary](dgx_pressure_overhead_failure_21869.json) binds all18
verification/task/log/scheduler records. The [capture receipt](dgx_pressure_overhead_capture_21869.json)
records complete terminal collection. Raw data remain under:

- `benchmarks/work/pressure_overhead_failure_21869/`.
- `benchmarks/work/pressure_overhead_scheduler_21869/`.

The historical submitted script hash remains recorded in the submission
receipt and git history; the corrected current script is not presented as
what ran. A regression check requires the pinned environment interpreter.
84focused launcher/plan/capture tests pass and shell syntax validation passes.
These checks do not claim successful inference or an overhead result.

Do not resubmit the existing plan against its occupied output directories.
A replacement must use a fresh complete18task plan/recipe/authorization,
retain native commands and budgets, and verify actual enumeration as well as
task selection before release. Preserve21869as a failed deployment. No
scientific setting or acceptance criterion is changed by this interpreter
repair, and no scientific timing has been admitted.
