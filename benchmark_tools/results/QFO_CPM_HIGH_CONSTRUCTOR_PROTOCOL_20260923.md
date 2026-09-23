# Failed High-CPM Constructor Diagnostic

## Scope Frozen Before Execution

The high-CPM replay 22081_1 failed by SIGSEGV in the fourth clustering worker
(index 3, profile_expanded). Its last marker precedes graph construction;
no optimizer-entry marker exists. The
[read-only audit](qfo_cpm_high_failed_payload_audit_20260923.json) checked
252 saved identities and all 25,501,180 edges over 984,137 vertices without
finding out-of-range endpoints, nonfinite/nonpositive weights or a changed
constructor-input digest. That evidence does not identify the native fault.

Run **one** new minimal-import Python-pair graph constructor on those exact
arrays. Use the existing pinned constructor diagnostic and endpoint/weight
comparison helpers. A fresh payload directory contains symlinks to read-only
inputs only; all observations are written to the new directory. The failed
payload and earlier completed stages remain untouched. The wrapper refuses
an existing output path, including a dangling symlink, and validates the
audit, helpers and saved records before and after its subprocess.

Enable Python fatal-error tracing from process start and retain the combined
worker log and native subprocess return code. The worker binds to one CPU,
constructs the undirected graph with explicit integer pairs, checks native
endpoints before weight assignment, then checks all weighted graph contents
against saved arrays. Neither OrthoHMM nor Leiden optimizer modules are
imported in this minimal worker. No partition, prediction or accuracy score
is produced. Full replay remains unauthorized regardless of this outcome.

## Resource and Interpretation Limits

Submit through `qfo_cpm_high_constructor_20260923.sh` with a clean frozen
executor commit: 1 CPU, 64 GiB, 2-hour limit, no requeue. GNU time is retained
as a diagnostic observation, not a controlled comparative timing result.
No DGX activity is involved. The original worker used a different import and
allocation history; this diagnostic is not an exact reproduction of those
conditions. There are no automatic repeats or alternate-format fallbacks.

Success means only that one fresh constructor preserved this graph. It cannot
rule out intermittent corruption, diagnose the historical fault, establish
optimizer correctness or authorize partial scientific output reuse. Failure
must retain the traceback/last observation and remain missing, not be hidden
behind retry. Review the result before designing any further intervention.

Fifteen focused tests pass, including a small fresh-process native graph,
endpoint/weight mismatch rejection, wrong-audit rejection, no-overwrite
protection and checks that the minimal worker did not import optimizer/core
modules. These do not establish behavior at the full failed-graph scale.
