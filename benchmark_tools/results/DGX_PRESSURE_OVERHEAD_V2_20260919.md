# Replacement Deployment After Pre-Inference Failure

The complete21869panel failed before native preparation because its system
interpreter lacked NumPy. Preserve that failure; it supplies no native
overhead measurements. This replacement uses all18tasks, not selected arms.

The new [v2 plan](dgx_pressure_overhead_plan_v2_20260919.json) has SHA-256
`b644e165dbf4d0beabf1cf4d9b6c314de522e3ebd1b91598ebebea99094c8fff`.
It is mechanically derived from the original pressure plan by relocating
`pressure_frontier_overhead_v1` to `pressure_frontier_overhead_v2` throughout,
then recording deployment_revision=2, supersedes_failed_array=21869,
the parent plan's path/hash, and the exact launcher interpreter:
`/home/jlsteenwyk/projects/orthohmm-publication/envs/orthohmm/bin/python`.
All prior plan fields otherwise match after this relocation, including
native commands, input identity/order, repetitions, observation mode,
resource limits, timeouts and numerical budgets. Tests verify this equality.

The launcher accepts only the explicitly pinned plans and now rejects an
interpreter differing from the v2 plan before entering the measurement
wrapper. No runtime package, frozen core, scientific setting or admission
criterion changed. The original recipe and occupied output root are retained.

The [pressure overhead protocol](DGX_PRESSURE_OVERHEAD_PROTOCOL_20260919.md)
applies unchanged. Before release, verify the fresh complete recipe, runtime
manifests, original inputs and actual frozen enumeration using the exact
launcher environment. Selection-only preflight was insufficient for21869.
Submit held; confirm durable controller capture and first poll; allow at
least60seconds after preparation; prohibit DGX reads during execution.

104focused tests pass, including all18v2selections, unchanged native work,
wrong-interpreter rejection, existing plan behavior and scheduler capture.
Execution still requires a separate exact recipe-bound authorization.
