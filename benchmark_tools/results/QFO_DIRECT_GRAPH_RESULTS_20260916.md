# Direct Graph-Stage Diagnostic

Job21326 completed0:0 in13:27, using frozen executorab364e4. Six fresh
one-CPU workers alternated minimal imports and frozen-worker imports. All used
the same saved graph and original int32 constructor input. No scientific worker,
optimizer, partition inference or accuracy scoring was executed.

| Import Mode | Repeat | Mismatches Before Weights | Mismatches After Weights |
| --- | ---: | ---: | ---: |
| Minimal | 0 | 6 | 6 |
| Frozen worker | 0 | 0 | 0 |
| Minimal | 1 | 0 | 0 |
| Frozen worker | 1 | 6 | 6 |
| Minimal | 2 | 6 | 6 |
| Frozen worker | 2 | 6 | 6 |

All constructor arrays matched the saved endpoints at both observations. Every
mismatching worker had the same six indices and replacements documented in the
[earlier construction diagnostic](QFO_CONSTRUCTION_DIAGNOSTIC_20260916.md).
The recorded complete difference reports were unchanged by weight assignment.

Independent admission checked terminal success, exact six-worker identity/order,
frozen executor/runtime, import isolation, common execution context, same-mode
module/library identity,278 file records and native-file/report agreement. It
reconstructed every post-weight native endpoint hash from the saved graph plus
the complete bounded mismatch witnesses. The pre-weight endpoint hash is only
implied by those witnesses: the executor did not separately record a full native
hash at that stage. This distinction is retained in the machine-readable report.

Verified snapshot: qfo_direct_graph_verified_20260916.json, SHA256
08352ac4af718e5b3a9e59d1cec2d81c3a333785aeeab8f22c1ed5975cfe7ff3.

## Interpretation

The mismatch is observed before weight assignment and without OrthoHMM or
Leiden imports in minimal workers. Thus neither weight assignment nor those
imports is necessary for this recorded failure. Both modes include clean and
mismatching workers; there is no evidence here that importing the scientific
worker fixes or causes the discrepancy. This does not identify a native-library
defect, hardware cause, or the cause of every historical partition difference.

The audit checks preserved observations, not the historical live graph object.
Direct setup and the extra observation alter memory allocation and timing.
The separate NumPy-versus-Python-pair diagnostic21327 remains running. No new
default or historical output substitution is justified; require native graph
integrity before using an optimizer replay as accuracy evidence.
