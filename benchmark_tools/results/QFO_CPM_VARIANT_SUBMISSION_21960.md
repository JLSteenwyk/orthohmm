# QfO CPM Variant Replay Submission

Submitted September 19, 2026 at 09:37:39 local time. Scheduler confirmed
array 21960, tasks 0-1, throttle 1, PENDING Dependency afterok:21958,
32 CPUs/192 GiB/bizon per task, 24-hour limit, Requeue=0.

Task 0 is cpm_low (0.08); task 1 is cpm_high (0.12). Both are fixed by the
September 19 parameter-neighborhood protocol. This is not an adaptive
parameter search. Inputs, HMM settings, seed and scientific runtime are
unchanged from the baseline. All four clustering/profile stages rerun.

Frozen executor: `benchmarks/work/publication_qfo_cpm_variant_v1`, commit
`2913d04d5289eeeee9d8f9940687494e96e7bb85` (pushed to origin/main).

Submission from the repository root:

```bash
sbatch --parsable \
  benchmarks/work/publication_qfo_cpm_variant_v1/benchmark_tools/results/qfo_cpm_variant_batch_20260919.sh \
  /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/publication_qfo_cpm_variant_v1 \
  2913d04d5289eeeee9d8f9940687494e96e7bb85
```

SHA-256 identities:

- `run_qfo_cpm_variant.py`:
  `60d29e3d8a0d97b95c9b141dabe7dc1125b5e98cd4a4a8d20e36645a49dbc548`
- Shared `run_qfo_cpm_control.py` worker:
  `71439aec3a8a22556b884242a18d066047712fe17a40da7300bf9fc2395cc102`
- `qfo_cpm_variant_batch_20260919.sh`:
  `7b05fd9dafa2fbc93350a238fe65044a80083a50baf52fddb438df6c867e785d`

Outputs: `benchmarks/results/qfo_parameter_cpm_replay_v1/cpm_low` and
`cpm_high`. Logs: `benchmarks/work/qfo_cpm_variant_21960_0.log` and `_1.log`.

Each run first requires completed control admission job 21958 and rechecks
its frozen source, context, authorization and file records. It then reruns
the independent validator from its frozen 40806f5 executor and requires
exact equality with the original control admission report before inference.
Failures persist without retry or overwrite. A failed upstream dependency
must be reported, not bypassed. Existing control/admission worktrees remain
unchanged by the shared worker's new explicit arm argument.

Validation: 253 focused tests passed in 12.76 seconds; batch syntax and
staged whitespace checks passed. The scientific import boundary was tested
in a real subprocess. No changed-arm output is yet admitted or scored.
Separate replay admission, candidates, phylogeny, conversion and benchmark
assessment are still required. Cached shared-host times are not controlled
efficiency evidence. Existing 11 Dependabot alerts remain unresolved.
