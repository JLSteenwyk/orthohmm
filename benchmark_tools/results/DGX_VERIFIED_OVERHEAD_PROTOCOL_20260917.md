# Prospective Verified-Wrapper Overhead Panel

The original six-run panel21640 is retained, including its passing numerical
budget and incomplete interpreter-environment identity evidence. This separate
repeat closes that identity gap prospectively and tests the GNU-time boundary.
No outcomes from this new panel have been inspected at protocol freeze.
The shorter engineering smoke21646 completed successfully and is not a paired
overhead result or scientific timing run.

## Fixed Design

Six fresh, sequential exclusive Slurm tasks on spark-7ff0, partition spark,
20CPUs96GiB, concurrency1,12-minute task limit,600-second native-command timeout.
CPU only. Indices0..5: sparse,sampled,sampled,sparse,sparse,sampled. Compare
sampled1/sparse0, sampled2/sparse3, sampled5/sparse4. No partial-outcome inspection,
optional stopping, selective retries or exclusion of unfavorable outcomes.
Preserve failures; any revised panel requires a separately recorded protocol.

Use unchanged native load fixture SHA256
80f97405b9361fa021b3857298ee7c4a3ee69e2f99d774cfba55a37afb0851ec,
8billion logical iterations per worker,20workers. Require all six native JSON
work/checksum records identical. Collector files remain the pinned
collector_load_recipe_v2 versions used in21640. GNU time wraps the native load
in both modes, with elapsed/user/system/max-process-RSS/exit fields. Require
zero native/GNU-time/collector exits and successful scheduler completion.

Sparse resource interval86400seconds, sampled1second; host interval30seconds
in both, retaining pre/post observations. Primary overhead statistic is
collector command-wall ratio minus1 for each fixed pair. Engineering budgets
remain median<=5% and every pair<=10%. Also require each native measured duration
>=60seconds, mean observed CPU>=18cores, at least3host scans for sampled and2
for sparse, no host observation errors, maximum observed foreign CPU<0.25cores,
and independently replayable complete raw collector records. These finite
observations are not a guarantee against missed short-lived/non-CPU contention.

## Identity and Boundaries

Verify runtime trees SHA256
2f38fc57683e51a7b6768b4709db16293590983640c41cfe12a56ed6036750ae
and system trees SHA256
4083d15c0ffc40013756c4588763aa8af24ca39165622665ee5074662e4b46ce
before and after each command. Additionally freeze the complete new recipe
directory (including wrapper, inventory helper, panel script and this protocol)
before submission; supply that manifest digest explicitly in the job arguments.
Its manifest lives outside the inventoried recipe directory to avoid circular
self-hashing. Record its exact digest in the progress ledger before submission.
The wrapper rejects changed manifests/files and records failed verification.

Identity checks execute in the collector's process, outside its measured native
command span. No extra shell remains in the task subtree. Record verification
wall time separately. Full-file hashing warms caches; neither arm is a
cold-cache experiment. The cgroup memory peak may include preparation. GNU-time
max-process RSS is not simultaneous process-tree RSS or cgroup memory. No
overhead subtraction is allowed.

Use an absent per-task Python cache prefix and disable bytecode writes; clear
PYTHONPATH/GOMP_CPU_AFFINITY/LD_PRELOAD/LD_LIBRARY_PATH/LD_AUDIT. Require absence
of /etc/ld.so.preload. Retain OMP_NUM_THREADS20,OMP_DYNAMICFALSE,OMP_PROC_BINDTRUE,
OMP_PLACEScores. Snapshot records exclude Python bytecode and git metadata.
System directory inventories capture file contents and resolved file symlink
targets, but are not a hermetic OS snapshot, nor proof against temporary changes.

## Interpretation

This is incremental periodic-observation cost for a small compute-only workload,
not a general overhead estimate or scientific accuracy/efficiency result.
Record GNU-time values alongside, not as a replacement for the prespecified
collector-wall statistic. Passing this panel does not itself authorize the27
scientific runs: the scientific preparation, per-run input checks, native-output
admission and failure inventory must still be integrated and tested.
