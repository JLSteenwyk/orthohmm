# Pending Workflow Retirement

## Scope

Live `squeue -u bizon` inspection on 2026-09-27 found 18 pending entries
and no running entries. These were old publication-workflow branches, not
evidence that additional analyses were progressing. All were cancelled with
`scancel --state=PENDING`, listing downstream jobs before prerequisites.
No running job, unrelated job, raw output, log or frozen executor was removed.
The DGX was not contacted.

## Branch Decisions

| Pending entries retired | Observed prerequisite state | Scientific disposition |
| --- | --- | --- |
| 21746, 21748, 21749, 21750, 21752, 21753, 21754, 21755 | 21746 user-held; downstream dependencies unfulfilled | Superseded original OrthoMCL branch. Retained results use the recovered-search pipeline and completed score admission 22176. Do not release the old 180-CPU inference job. |
| 22082_1, 22084_1, 22086_1, 22088_1, 22090_1, 22092_1, 22094_1, 22096_1 | 22081_1 FAILED, exit 1:0; first admission DependencyNeverSatisfied | High-CPM experiment remains failed/unavailable. Cancelling its unstarted descendants does not supply scores or repair inference. |
| 22105_14 | 22103_14 TIMEOUT after 1-10:51:21; DependencyNeverSatisfied | Old BLAST batch validator is superseded by replacement admission 22161, which completed 0:0. Original timeout remains a failure record. |
| 22156 | 22155 FAILED, exit 1:0; DependencyNeverSatisfied | Recovered high-CPM candidates remain unavailable. No retry or admission was launched. |

Retained OrthoMCL evidence is
[score admission 22176](qfo_recovered_score_admission_22176.json) and the
[eight-method manifest](qfo_corrected_comparison_20260926_v7/manifest.json).
The latter binds this admission and recovered pair conversion, not the old
21746 branch. Accounting confirms 22176 COMPLETED 0:0 in 00:06:02 and
22161 COMPLETED 0:0 in 00:00:14. Intermediate search admission 22163 FAILED
1:0 remains in accounting; the existence of later admitted results does not
rewrite that attempt as successful.

## Verification

Cancellation command, exit 0:

```sh
scancel --state=PENDING 21755 21754 21753 21752 21750 21749 21748 21746 22096_1 22094_1 22092_1 22090_1 22088_1 22086_1 22084_1 22082_1 22105_14 22156
```

Immediate `squeue -u bizon` returned only its header. `sacct -X` reported
17 entries CANCELLED by 1000, all with elapsed 00:00:00 and end time
2026-09-27T11:49:30 (scheduler-displayed time). The missing array entry
22086_1 was checked directly with `scontrol show job 22086_1`: CANCELLED,
RunTime=00:00:00, NodeList empty and AllocTRES=(null), same end time.
The cancellation-time StartTime shown for some array records is not evidence
of execution. Cancellation exit code 0:0 is not a successful scientific run.

This verifies an empty account queue at inspection time, not absence of
non-Slurm processes or other users' workloads. It does not establish a quiet
timing host. Dedicated timing, appropriate uncertainty, missing TreeFam source
data and publication packaging remain unfinished. Resuming a failed branch
requires a justified new execution plan, not release of these obsolete jobs.
