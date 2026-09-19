# DGX Read-Crossing Control Result

All three prespecified controls completed in job 22018, exit 0:0, elapsed
three seconds. Slurm reserved all 20 DGX CPUs exclusively with 4 GiB memory;
the observer step was bound to CPUs 0 and 1. No native timing panel was live.
Every owned finite service used CPU 0, finished, and had its cgroup removed
between the crossing snapshot's root read and next ancestor read.

| Trial | Service process CPU (s) | Callback wall (s) | Full-span root-minus-observer (CPU s) | First partial span (CPU s) | Second partial span (CPU s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.750202 | 0.818291 | 0.809154 | -0.010846 | 0.820000 |
| 1 | 0.750191 | 0.827875 | 0.814570 | -0.011430 | 0.826000 |
| 2 | 0.750238 | 0.829050 | 0.811985 | -0.011016 | 0.823001 |

All fixed full-span and nested lifecycle response checks passed. The raw
comparisons, service command/limits, affinity, membership, observed deletion,
scope/boot identity and callback ordering reproduce independently. Negative
partial-span differences are retained: a root counter sampled before an
event and other counters sampled later do not yield atomic differences.
This is direct evidence of read-window skew in a known engineering control,
not a causal explanation of native-run residuals.

## Provenance And Replay

Deployed ten files exported from committed source `6599c6e`. Local and DGX
archive SHA-256 matched:
`a0f3c2332b0971229273a148b4d902b7c9755391120b0508ea36cddcb358ccb9`.
Byte comparison of every file and the complete source inventory passed;
the output directory was absent before submission. Tar's metadata comparison
reported expected UID/GID differences from extraction as the remote user,
not content differences. No prior control directory was overwritten.

Submitted the committed `run_dgx_lineage_read_crossing.sh` with loader and
Python-path environment overrides removed and TMPDIR=/tmp. The captured
terminal scheduler record and accounting are retained beside 35 raw files
(160,395 bytes) in `lineage_read_crossing_22018/`. Kernel identity was
6.11.0-1016-nvidia, aarch64, with Python 3.12.3. Source hashes were unchanged
across execution and independently matched the frozen Git objects.

Retained replay: `lineage_read_crossing_replay_22018_20260919.json`, SHA-256
`cbd3a1aca685477554b951742523057756b1e06d6631cc0484c4abc508e6bd89`.
Original control report SHA-256:
`582717af59d3a8713c8a3c7cafa2a43b54b79fced7a808a790c338622a124641`.

```sh
python -B -m benchmark_tools.replay_lineage_read_crossing_control \
  --directory benchmark_tools/results/lineage_read_crossing_22018 \
  --scheduler benchmark_tools/results/lineage_read_crossing_scheduler_22018.txt \
  --repo . --output /tmp/lineage_read_crossing_replay_22018.json
```

Use the retained checkout/source identity and an absent output path. The
replay validates all three trials, not only successful selected records.
All 62 focused tests pass, including real-trial replay and rejection of
changed commands, service failures, missing deletion, wrong scopes, event
ordering, altered results, missing files, changed sources and resource
allocation. No scientific executor or running local analysis was changed.

## Scope

This validates bounded aggregate observation during one deliberate owned
user-service lifecycle per trial. It does not validate all forms of churn,
individual-process attribution, false-positive rates, root/user-slice
specificity, callback-enabled overhead or absence of non-CPU interference.
The injected read delay is not used in scientific measurements. No original
flag is removed, no time corrected and no historical run promoted.

Next: prospectively separate root/user-slice and native-only accounting
under known workloads, then define scientific inclusion before controlled
resource comparisons. The full publication goal remains incomplete.
