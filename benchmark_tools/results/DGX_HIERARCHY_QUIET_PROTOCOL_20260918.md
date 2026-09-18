# Native Hierarchy Quiet Control

One prospective complete three-method control following the known operator
SSH/SCP overlap in21817. Retain that array, its flags and all previous runs.
This is not selective repetition until a threshold passes. Run all three
methods once in the same fixed order, even if earlier tasks fail or flag.

Reuse the frozen645-protein/eight-species specification, runtime/system
manifests, native commands, hierarchy collector and0.25core/negative0.5CPU-
second/steal/gap checks. Only relocate outputs to hierarchy_quiet_smoke_v1
and prepare a distinct hierarchy_quiet_recipe_v1. Use sequential exclusive
spark-7ff0 tasks,20CPUs/96GiB,900second native timeout and no requeue.

Complete all DGX file transfers and recipe identity checks before submission.
Submit with a60second future eligibility time to let the submission SSH
session exit before any native task starts. While the array is pending or
running, issue no DGX SSH, SCP or other remote commands. Poll only the shared
scheduler from the main host. Do not archive, inspect remote logs, change
services or stop unrelated workloads until all tasks are terminal. Retain
the scheduler BeginTime/StartTime and the submission-session completion
observation. This operator policy does not prove other activity absent.

After terminal completion, archive all outputs and repeat unchanged native
and counter validation. Report every flag and signed hierarchy residual.
Compare descriptive patterns to21817 only; do not estimate causal SSH
overhead or tool speed from two small, ordered arrays. A clean result would
support only this quiet-control observation, not general isolation, native
overhead calibration or scientific timing admission. A flagged result is not
grounds to raise the threshold or selectively rerun a method.
