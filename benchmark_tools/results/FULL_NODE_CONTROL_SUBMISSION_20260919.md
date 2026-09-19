# Full-Node Control Submission

## Pinned Deployment

Exported committed sources from `dd912c6`, not the dirty worktree, using
`git archive`. The447 regular files in the deployed recipe were compared
byte-for-byte with their Git blobs through SHA-256 and size checks. No
external symlinks were present.

- Archive: `benchmarks/work/full_node_controls_recipe_dd912c6.tar`.
- Archive SHA-256: `d2d91410ebb454804edf6084914d4798f27f2e9ca4b801f6f16bc6b915dda952`.
- Remote recipe: `/home/jlsteenwyk/projects/orthohmm-publication/full_node_controls_recipe_v1`.
- Manifest: `full_node_controls_recipe_20260919.json`.
- Manifest SHA-256: `796da266cf839441d239b1ffae4c372ff616b8feacd64dada84746f0d6c31e62`.
- Output: `/home/jlsteenwyk/projects/orthohmm-publication/full_node_controls_v1`.
- Protocol SHA-256: `76ebeb56b1b7874aa03528064d77fd1c38ec1151a74095204a8393632b0e834d`.

Both existing application and system runtime inventories verified on the
DGX before submission; the driver verifies their pinned hashes and contents
again before and after the panel. The copied batch file matches the committed
wrapper exactly. A controller queue check found no DGX jobs; the two current
one-second `vmstat` observations reported100% idle, zero swap activity and
zero steal. This short preflight does not prove isolation during execution.

## Failed Launch Preserved

Job21916 failed at shell startup with exit1:0, zero scheduler elapsed time,
before Python or any control trial. The batch inherited the local
`LD_LIBRARY_PATH` CUDA override, which its explicit guard rejects, plus local
working/TMPDIR paths unavailable on the DGX. The log records the path
warnings; the scheduler record and local exported environment support this
diagnosis. No panel output directory existed after failure.

Recorder21917 completed successfully and retained its terminal record.
Evidence: `full_node_control_capture_21917_20260919.json`,
`full_node_control_scheduler_21916_20260919.txt`,
`full_node_controls_21916.log`. This is a failed submission, not a failed
control outcome. It is not omitted or overwritten.

## Corrected Submission

Submitted the same batch and recipe as21918 with explicit launch options:

```text
--chdir=/home/jlsteenwyk/projects/orthohmm-publication/full_node_controls_recipe_v1
--export=ALL,LD_LIBRARY_PATH=,LD_PRELOAD=,LD_AUDIT=,PYTHONPATH=,TMPDIR=/tmp
```

No workload duration, condition order, threshold or protocol changed.
Controller confirmed21918 RUNNING at00:00:11 on spark-7ff0, exclusive20CPU,
96GiB,45-minute limit, zero restarts and no requeue. No SSH inspection is
performed during this live panel. Native commands retain the60-second timeout.

Recorder21919 runs separately on bizon with1CPU/512MiB and a one-hour limit.
It uses an extracted copy of the same committed recipe, recording job21918
under `benchmarks/work/full_node_control_scheduler_21918/`; its log is
`benchmarks/work/full_node_control_capture_21919.log`.

After terminal capture, collect the archive and validate runtime records,
scheduler resources, all nine trials, raw measurements and worker witnesses.
Report failed controls and missed detections without selective retries.
These engineering controls do not themselves admit the27 scientific scaling
runs, establish monitor overhead, identify historical foreign workloads or
prove publication readiness.
