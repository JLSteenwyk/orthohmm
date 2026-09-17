# QfO Construction Diagnostic

Job 21323 completed successfully in 18:44. This construction-only diagnostic
never invoked the Leiden optimizer and produced no accuracy scores or partitions.
Frozen executor: f2827a6. It alternated three original-int32 and three explicit
int64 construction workers using the same saved graph and one-CPU setting.

| Worker | Constructor Input | Native Endpoint Mismatches |
| --- | --- | ---: |
| 0 | Original int32 | 0 |
| 0 | Explicit int64 copy | 0 |
| 1 | Original int32 | 0 |
| 1 | Explicit int64 copy | 6 |
| 2 | Original int32 | 0 |
| 2 | Explicit int64 copy | 0 |

All original and converted constructor arrays matched the saved endpoints. The
mismatching native graph differed at edge indices 23493880 through 23493885.
The original expected pairs linked vertex 828439 to vertices 851595, 855360,
885964, 922655, 924806 and 931685. Observed replacement endpoints were 65163,
0, 0, 0, 0 and 917504. Native edge tuple and source/target access agreed, while
get_eid did not find the six expected pairs. Vertex count, edge count and the
complete ordered weight fingerprint were unchanged.

Independent report admission checked terminal scheduler success, frozen executor
and runtime, all six worker identities and saved inputs, 257 file records, and
agreement between report entries and preserved native observations. Applying the
complete six-edge witness list to the saved endpoint stream reproduces the full
reported native endpoint SHA256. All five clean workers' hashes reconstruct too.
Snapshot: qfo_construction_verified_20260916.json, SHA256
499f90e06b6da311356c80f435149ba0cf132f75b13e2eacf720e02ea8f2382e.

## Interpretation And Next Steps

Explicit int64 conversion alone is not sufficient to eliminate this observed
failure. This is not evidence that int32 is generally correct, that the conversion
caused the failure, or that a particular native library or hardware component is
responsible. The audit validates preserved observations and hash consistency;
it does not independently inspect the historical live graph object.

Next isolate direct graph construction from the frozen worker's surrounding
imports and setup, and inspect graph identity before and after weight assignment.
Keep optimizer execution separate until native endpoint integrity can be enforced.
Preserve failed and clean observations; do not select a partition by accuracy or
closeness to historical output. QfO replay-based accuracy comparisons remain
unresolved pending a reproducible, correctly constructed baseline.
