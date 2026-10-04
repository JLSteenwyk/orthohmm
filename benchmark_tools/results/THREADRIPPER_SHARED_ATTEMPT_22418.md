# Shared-Host Attempt 22418

After the independently reviewed index-21 checkpoint is committed/pushed at
`67727254`, the existing prepared continuation validates all 22 earlier
identities and releases only index **22** as **22418**: phylogenetic OrthoHMM
`satellite_v2`, eight proteomes, repeat 2. No failed identity is retried,
scientific configuration changed or unrelated workload disturbed.

The [public submission](threadripper_shared_submission_22418.json) matches the
canonical launch receipt: 1,179 bytes, SHA-256
`1e110b3d9edc4b02c056141c1e468c90f3363e80322b58a997a419d11ad74b19`.
Request: 7,765 bytes, SHA-256
`8e6f57d09a26712ee819478f67879f638778d7bf56cd915180042900058ac989`.
The retained preparation/source/resource/policy/readiness bindings are reused;
no completed calibration, regression or runtime refresh is repeated.

## Actual Native Handoff

The one-shot [live snapshot](threadripper_shared_live_22418.json) passes:
16,770 bytes, SHA-256
`01b4724a6c8bf6b860c08a9d065176ed53eff9d4c8fa13ab1b10a02227bb6198`.
It validates the running native PID and job-step membership, request digest
in the scheduler comment, passed preflight, true go gate, environmental release
and prepared-before-request ordering. Preparation precedes release request by
59.203212230s; request-to-bound-review is 9.816271440s. Native affinity remains
physical CPU IDs 0-31, under a shared 64-SMT-slot Slurm allocation and enforced
137,438,953,472-byte memory limit. The allocation retains 26h and no requeue.

Both preflight available-memory observations, 754,510,417,920 and
753,059,823,616 bytes, exceed the unchanged 128-GiB safety floor. Recorded
foreign CPU demand is 39.2836691948 core equivalents. The first observed
process-start gap is 30.000301076s, within the retained cadence bound, but
does not establish whole-run coverage or causal observer overhead.

The actual frozen command includes reconciliation, inferred species tree,
MAFFT, FastTree, satellite_v2 candidates, species-overlap gene-tree rooting,
positive-paralogy pairs and minimum-variance species-tree rooting. The startup
log currently shows its built-in HMM/k-mer search on eight proteomes. This
establishes native execution, not completed phylogenetic inference or accuracy.

## Continuation Boundary

Fresh accounting confirms parent/batch/native step RUNNING. The launcher and
snapshot tool handles are terminal; the job is not restarted. Existing
[table v22](threadripper_shared_panel_snapshot_20261004_v22/panel.json) and
figure v21 remain unchanged: 22 reviewed, 20 measured, 19 eligible, exclusions
0/17/20 and two complete three-eligible-repeat cells. Five identities remain
unreviewed: this one live and four not submitted. No new resource endpoint,
success summary or accuracy score is admitted from startup evidence.

On resumption poll this exact job, retain terminal controller before purge
when available, and use the existing prepared terminal reviewer only after
actual termination. Do not rerun this one-shot snapshot or successful
preparation/diagnostics without a concrete reason. Launch index 23 only after
independent review advances history. Preserve all attempts and failures.
Timings remain shared-host observations with unknown, potentially
method-dependent contention distortion, not isolated efficiency evidence.
Final panel, manuscript/resource reconciliation and versioned/archive release
remain incomplete; the full publication goal stays active.
