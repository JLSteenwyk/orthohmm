# Checked Replay Environment Failure

Job21329 terminated FAILED1:0 after46m36s. The full replay is not admitted and
must not be scored or substituted for the historical publication results.

Initial and multipass payload checks completed with976,504genes and respectively
349,898 and308,190groups. The third clustering call, profile_base, returned from
the optimizer but failed the subsequent environment gate: OMP_NUM_THREADS was32
instead of1. CPU affinity remained a single CPU (8). OPENBLAS_NUM_THREADS and
MKL_NUM_THREADS remained1. Frozen profile_expansion.py:702 assigns the requested
profile CPU count to os.environ, which subsequent subprocesses inherit.

All three saved native-boundary observations report equal before/saved/after
graph fingerprints and optimizer_returned. This does not waive the environment
failure or establish independent full-replay admission. The fourth clustering
call was never reached. The earlier NumPy-constructor corruption investigation
remains a separate issue; this failure is not evidence that it recurred.

Preserved snapshots:

- qfo_checked_full_replay_failed_20260917.json, SHA256
  6ed8d21b2d415c8b80c1069f205ea461bc290462d5ea4d24e9fdebb17647aa94.
- qfo_checked_full_replay_failed_worker_20260917.json, SHA256
  b20e33e51e01334084b1a95729ee1dd11dcfe31b22c0acc60bb930490ad381cb.
- Raw logs, payloads and partial outputs remain under
  benchmarks/results/qfo_checked_full_replay_v1 and its job21329 scheduler log.

Correction: the executor supplies explicit OMP/OPENBLAS/MKL thread values of1
to each clustering child. It records the inherited values and overrides. The
parent environment, profile32CPU settings, scientific sources, arrays, optimizer
settings, affinity, and validation criteria are unchanged. A regression test
covers all four calls with inherited32/8/4 thread settings and proves the parent
environment is not mutated. Environment errors now list the offending keys.

A separately named v2 replay uses fresh outputs. The full run is repeated
because the frozen replay has no validated resume point at the failed callback;
reusing partially consumed mutable state would introduce another unvalidated
execution route. The old outputs are retained for stage-by-stage comparison.
This retry is driven by a diagnosed execution defect, not accuracy selection.
The current independent admission script pins v1/job21329/executor96333fd and
must be explicitly updated for v2 provenance before any completed v2 result is
admitted. No native scientific library or frozen core is patched.
