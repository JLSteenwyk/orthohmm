# Bounded Session Follow-up For Root Context Controls

## Rationale And Scope

Panel 22019 remains an incomplete experiment: three validated native-only
conditions, one failed owned-service launch and eight unrun conditions.
Its archived journal shows the user manager shutting down ten seconds after
the submitting SSH connection closed. No original condition is replaced.

Run one separately identified complete 12-trial panel with the same order,
workloads, timeouts, resource limits, counters and response criteria frozen
in `ROOT_CPU_CONTEXT_CONTROL_PROTOCOL_20260919.md`. The only intended
execution change is maintaining the submitting login session through job
termination. Use fresh `root_context_controls_recipe_v2` and
`root_context_controls_v2` directories and a separately pinned source recipe.
Do not pool missing repetitions from 22019 with this follow-up.

## Bounded Session

`submit_root_context_session.py` checks and retains an empty DGX queue, then
opens one noninteractive SSH connection running `sbatch --wait --parsable`.
[Slurm documents](https://slurm.schedmd.com/sbatch.html#OPT_wait) that `--wait`
does not exit until the submitted job terminates and returns its exit code.
The remote client is bounded by GNU timeout at 1,020 seconds plus a ten-second
termination grace, and the local observation by 1,050 seconds. The scheduled
panel retains its 900-second allocation limit. A timeout or connection loss
is an observation failure: inspect the existing Slurm job, do not assume it
ended, cancel unrelated work, or resubmit.

The held session carries no interactive commands, file transfers or host
inspection during the panel. The remote sbatch client waits on Slurm; it is
not an inference worker. Retain this additional session/client presence as
part of the control environment rather than claiming zero observer work.
No persistent lingering setting or unrelated user service is changed.

A bounded preparatory SSH probe retained at
`root_context_held_session_probe_20260919.json` observed the same user-manager
scope before and after a 20-second wait. This is evidence for the session
lifetime mechanism only, not a complete-panel or overhead validation.

## Evidence And Failure Handling

Retain local queue, launch and completion/timeout receipts, the SSH command,
its exit status and unambiguous job identifier when available. Independently
retain terminal scheduler evidence and retrieve the entire panel only after
termination. Verify source/runtime identities and all workload witnesses
using the v2 audit binding. Check historical journal evidence that the user
manager did not stop during the job. A user-manager lifecycle failure or
unconfirmed service cleanup invalidates affected conditions and stops further
work as before. All failed and unrun outcomes remain explicit.

Repeat the prescribed all/common interval descriptions and fixed-block
differences without treating intervals as independent replicates. Require
the unchanged five-CPU-second aggregate response in each valid user-contended
condition, retaining low responses. This remains an observer-engineering
experiment, not native-tool performance evidence. Native overhead, causal
attribution, non-CPU isolation and scientific timing admission remain separate.
