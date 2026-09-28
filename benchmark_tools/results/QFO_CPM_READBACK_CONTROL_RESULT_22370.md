# Isolated Readback: Crash Not Reproduced

The single control prescribed in
[the protocol](QFO_CPM_READBACK_CONTROL_PROTOCOL_20260928.md) completed as
Slurm job 22370, exit 0:0, seven seconds, one CPU and 8 GiB. Protocol, runner
and ten passing tests were committed and pushed as e4e84af5 before launch.
There was no retry or scientific inference. This is not controlled timing.

The retained Python executable ran with site initialization disabled,
`PYTHONMALLOC=debug`, faulthandler, and default enabled GC thresholds
700/10/10. The original reader function was extracted unchanged from the
frozen module AST; module-level scientific imports were omitted. It read
the initial partition and the failed diagnostic's saved refined partition.
All 984,137 genes appeared exactly once, in 390,845 refined groups.

| GC generation | Collections during the two readbacks |
| --- | ---: |
| 0 | 916 |
| 1 | 83 |
| 2 | 7 |

No named scientific package was loaded; the receipt retains the full module
list, not only that check. Child stderr was empty. Eleven input/source/config/
output records were independently rehashed after terminal completion, and
runner bytes match the committed source. The [receipt](qfo_cpm_readback_control_22370.json)
includes scheduler/controller evidence and raw-artifact identities.

This control shows successful parsing and garbage collections for the actual
saved partition bytes in an isolated process. It does **not** reproduce the
heap state or allocation history after refinement, establish memory safety,
implicate a particular extension, rule out intermittent interpreter/hardware
faults, or show equivalence of all runtime libraries. No root cause or fix
has been identified. Admission 22155 remains failed; high-CPM accuracy stays
missing with the original multiplicity handling intact.

Further debugging should isolate the preceding import/data-loading/refinement
state or obtain native evidence from a reproducing failure. Repeating this
successful isolated reader, repeating the startup control, or substituting
this result for full scientific admission would not resolve the requirement.
Any stage-isolation experiment needs a separate bounded protocol.
