# Dual-Bracket CPU Control Result

## Execution and Validation

Job 21911 completed on spark-7ff0 with exit 0:0 in 19 seconds, without
requeue or restart. The [scheduler record](dual_bracket_controls_21911_scheduler.txt)
shows two requested task CPUs, 256 MiB and exclusive-node allocation;
Slurm consequently allocated all 20 node CPUs. The worker and observer used
distinct permitted CPUs, with the competitor pinned to the worker CPU.
This is not a two-CPU concurrent shared-node timing result.

The node was idle with CPUAlloc=0 before submission. No unrelated process
or service was stopped. Recipe source commit: `081437b`. Archive
`benchmarks/work/dual_bracket_recipe_v1.tar` SHA-256:
`a13c2a2ca14c884b3f1865e3808b64ca8a35a2006b8afffa1e69ac500c45f9ea`.
The transferred archive hash matched, and all 433 extracted files matched
the archive before dispatch. Python was 3.10.13 (conda-forge).
The protocol was frozen in the preceding commit, before these outcomes.

After terminal scheduler state, copied the complete trial directory to
`benchmarks/work/dual_bracket_controls_21911/`. Independently verified all
432 reported source/protocol hashes against both the transferred recipe
archive and local sources. For every trial, compared standalone before,
after, ready, done and trial JSON files with the assembled report; checked
competitor exit and parsed stdout where present. Replayed both CPU brackets,
native pressure and workload validation, then reproduced the entire summary
exactly. All nine trials are valid and retained.

## Prespecified Results

| Trial | Mode | Narrow residual (cores) | Wider residual (cores) | Narrow expectation |
| --- | --- | ---: | ---: | --- |
| 0 | Quiet | 0.064467 | 0.075421 | Pass |
| 1 | Native-only | 0.097951 | 0.116202 | Pass |
| 2 | Contended | 0.514531 | 0.520019 | Pass: positive CPU flag |
| 3 | Native-only | 0.099655 | 0.118244 | Pass |
| 4 | Contended | 0.514224 | 0.519707 | Pass: positive CPU flag |
| 5 | Quiet | 0.108187 | 0.119146 | Pass |
| 6 | Contended | 0.514333 | 0.519819 | Pass: positive CPU flag |
| 7 | Quiet | 0.074293 | 0.079709 | Pass |
| 8 | Native-only | 0.118447 | 0.127743 | Pass |

All six uncontented narrow screens pass; all three injected-load screens
flag positive excess CPU. Native CPU PSI `some` increases relative to the
same block's native-only control by 702,387, 703,756 and 699,797 microseconds,
respectively, exceeding the prespecified 100,000-microsecond requirement.
Injection CPU dose, scope, affinity, overlap and post-work observation delay
all pass validation. No trial was excluded or repeated.

## Limits and Next Step

Both wider and narrower screens respond correctly in this single-core
experiment; it does not establish that the narrower collector is universally
better. The earlier archived replay identified wider-bracket inflation under
real native workloads, but 46 narrow residual flags remained. Their causes
and full-node native behavior still require investigation.

This result does not establish observer overhead, memory/I/O sensitivity,
an interference bound, or quiet full-node inference. Cgroup identity churn
still fails validation. No thresholds or historical outcomes changed, and
0/27 scientific scaling runs are admitted. Next evaluate the narrow bracket
with native/full-node workloads under a separately pinned recipe before
the complete overhead panel and scientific timing admission.

## Retained Evidence

- [Compressed report](dual_bracket_controls_21911.json.gz), SHA-256
  `b665855013c44a2f8da899934abaeecb5d742ece2d0a83e8a3ea1aa6a83e6685`.
- Decompressed report SHA-256:
  `806267722d3a77d6188de30805ff2f4b86bfe159d6ec41a7a4f362d32b2969c9`.
- Scheduler record SHA-256:
  `6dde87bbb28279a8bd286e395253e46a3632799deeb9554c7ec4112abd3b08b3`.
- Runner SHA-256:
  `d209e6a1bae11afb264f87c6a3ab12ae4f356558e756b9e45fa348b88dba6d60`.
- Batch script SHA-256:
  `89bf2f14596c2e4d41e927339d95895001f2004f849a1287da692232e73be4f9`.

Runner, reader and existing pressure controls pass 47 focused tests;
`bash -n` passes for the batch script. Reproduce report validation using
`run_dual_bracket_controls.validate_trial(row, report['job_id'])` for all
nine rows and compare `summarize(report['trials'])` to `report['summary']`.
These validation functions include raw dual-bracket and pressure replay.
