# QfO CPM Control Replay Submission

Submitted September 19, 2026 at 09:24:07 local time. Scheduler inspection
confirmed job 21956 PENDING (Priority), 32 CPUs, 192 GiB, bizon, 24-hour
limit, Requeue=0. Submission is not successful execution or admission.

Frozen executor: `benchmarks/work/publication_qfo_cpm_control_v1`, detached
commit `ef5fcd6f9d0ed46d8facc67631dcc7d3d3bb8f3e` (pushed to origin/main).
Root is `/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm`.

SHA-256 identities:

- `benchmark_tools/run_qfo_cpm_control.py`:
  `836da56448214cabc945dcdaec84d38dd09222189b080424964713db908b9002`
- `benchmark_tools/results/qfo_cpm_control_batch_20260919.sh`:
  `bc8dd05506a4db6f9751402fc6f4a3d6180680a992053e69027fe103388de968`

Submission from the root:

```bash
sbatch --parsable \
  benchmarks/work/publication_qfo_cpm_control_v1/benchmark_tools/results/qfo_cpm_control_batch_20260919.sh \
  /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/publication_qfo_cpm_control_v1 \
  ef5fcd6f9d0ed46d8facc67631dcc7d3d3bb8f3e
```

Log: `benchmarks/work/qfo_cpm_control_21956.log`.
Output: `benchmarks/results/qfo_parameter_cpm_replay_v1/control`.
The unchanged 0.1 control must reproduce complete memberships at all four
baseline stages. Worker payloads are checked independently during replay;
the final control still requires separate admission before changed CPM
arms can proceed. Neither accuracy nor controlled end-to-end timing is
established by this job. No other jobs were stopped or resource caps changed.

Pre-submission validation: 137 focused tests passed in 6.06 seconds, batch
syntax and staged whitespace checks passed. The remote reported 11 existing
Dependabot alerts (3 high, 7 moderate, 1 low); these remain unresolved.
