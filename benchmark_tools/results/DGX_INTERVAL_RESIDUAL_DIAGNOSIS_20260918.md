# Native Interval Residual Diagnosis

Retrospective diagnostic of all three retained21810native smokes. This does
not change the operational screen, rerun any method or upgrade timing evidence.
The original counter records are fixed by their SHA-256 identities; every
original interval and whole-command screen was replayed exactly before the
new arithmetic. See [machine-readable result](dgx_interval_residual_decomposition_20260918.json).

## Findings

The high-sensitivity smoke has no flagged interval. For the other two,
zero-based interval3 remains flagged:

| Configuration | Original host-minus-native CPU seconds | Outer-minus-inner host CPU seconds | Observer leaf CPU seconds | Inner-host-minus-native CPU seconds |
|---|---:|---:|---:|---:|
| satellite_v2 | 0.324697 | 0.010000 | 0.002877 | 0.314697 |
| OrthoFinder full | 0.269909 | 0.010000 | 0.003075 | 0.259909 |

The outer-minus-inner difference equals the summed counter increments within
the two endpoint read windows. It is an exact telescoping identity of the
reported host counters, not a physical bound on CPU usage or accounting lag.
The observed endpoint increments and observer leaf usage are small relative
to the flags. Neither observation explains away the discrepancy. The inner
residuals remain positive; all negative residuals elsewhere are also retained.

Native aggregate versus host user/system differences are reported for every
interval, but do not identify the source of the remaining CPU. The observer
leaf is only the batch user task, not the full batch step, job container,
Slurm daemons, asynchronous kernel work or unrelated system processes. Its
CPU usage must not be represented as total observer overhead.

## Implication For Next Measurements

Do not loosen0.25cores, drop flagged intervals, subtract residuals from wall
time or selectively repeat these methods. To distinguish unobserved job work
from outside-job activity, prospective controls should capture the native-step,
batch-step and job-parent aggregate counters with explicit hierarchy and read
brackets. Parent-minus-children discrepancies require checks for overlapping
scope and counter-accounting delay; they are not automatically foreign load.
Controls need sustained and completed-burst loads plus read-only/no-observer
comparisons before a scientific inclusion policy can be justified. Those
prospective controls are not implemented or executed by this diagnosis.

No claim of exclusive non-CPU resources, negligible overhead or controlled
comparative speed follows. The original27timings remain descriptive and the
scientific timing gate remains unmet.

## Verification

58focused tests passed across the new diagnostic and existing interval,
complete-command and bracketed screens. Thirteen new tests cover raw counter
parsing, negative/decreasing counters, exact read-window algebra, preserved
scope/error guards, independent tick arithmetic on retained quiet/burst
controls, mismatch rejection and retention of both adverse smoke results.
No frozen observer source was modified.
