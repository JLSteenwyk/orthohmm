# Prospective DGX Lineage Lifecycle Control

Run three sequential, finite engineering controls in an exclusive spark
allocation with a two-CPU-bound observer step and 4 GiB memory. The exclusive
allocation reserves all 20 node CPUs; two is the observer step's affinity,
not the node allocation size. This is not tool timing, collector
overhead admission or a replacement for the failed panel 21920.

The observer remains in its Slurm step cgroup. Snapshot the stable user-manager
ancestor and the root-to-observer lineage. Start a uniquely named owned user
service `orthohmm-lineage-JOB-INDEX.service` for each control, with one CPU
chosen from the allocation affinity, RuntimeMaxSec=10s, MemoryMax=256M and
TasksMax=8. Its only payload is the existing finite .75-process-CPU-second
burn function. The service is deliberately outside the Slurm job to test
outside-job accounting. It is not unrelated background activity.

Record argv, service stdout/stderr/return code, process CPU, PID, membership,
affinity and checksum. Require .75 <= process CPU < .85 seconds, exactly one
selected CPU, a service directly beneath the expected manager subtree,
and membership outside the observer scope. Wait at most 15 seconds for the
service cgroup to disappear, recording each existence check. Only then can
the second snapshots support the lifecycle test. Preserve a failure if it
does not disappear. No unrelated service may be stopped; these transient
units use `--collect` and a finite runtime cap.

For each of three controls, require a manager aggregate delta of at least
.5 CPU seconds after service deletion and a signed root-minus-observer delta
of at least .5 seconds. These fixed engineering response thresholds follow
the earlier .75-second workload controls; they are not comparative timing
eligibility thresholds or estimates of causal interference. Keep all signed
differences. Other system activity may contribute to either delta, so this
does not establish exact recovery of workload CPU or false-positive rates.

Retain all three outcomes without selective retries, all raw snapshots,
service lifecycle observations, source hashes, kernel/Python identity and
Slurm allocation/terminal records. Failed observations retain partial
evidence. Require stable boot/lineage identity and unchanged observer
membership and sources. Do not admit results based only on a process exit.

This checks completed descendant accounting at an existing user-manager
ancestor and an outside-observer response on this DGX kernel. It does not
test the exact system.slice service identities in the historical failures,
service churn during a read, attribution specificity, full-node inference,
non-CPU isolation or overhead. Native integration/replay, a complete new
overhead panel and a prospective scientific inclusion policy remain separate.

## Preflight Failure And Launch Correction

Attempt 21988 failed before output-directory creation or service launch:
the exclusive batch allocation exposed all 20 CPUs and the two-CPU affinity
guard rejected it. The corrected wrapper uses an explicit exact two-CPU
Slurm step with core binding, a remote-valid working directory and TMPDIR.
Sources/output paths are new v2 locations; v1 and its failure remain intact.
All three controls will run afresh with unchanged workload and response
thresholds. This is a launch correction, not selective outcome-based retry.
