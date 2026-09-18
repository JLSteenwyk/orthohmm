# Corrected QfO All-Hit Graph Execution

Job21813 completed0:0 in37:34 onbizon with32CPUs and384GiB allocated.
The [source-bound execution report](qfo_sequence_all_hits_execution_21813.json)
has SHA-256
`e91752ab1c2a12837c2c663eb55c2cb80b1b69db10ceb6f4e3f96946332042fb`.
Its status remains `sequence_checked_graph_complete_pending_admission`.

The frozen graph executor reports completed initial and multipass checked
clustering, unchanged runtime identities and both complete984137-gene
partitions. Reported groups:217059multipass and417273multipass_refined.
These are unscored output counts, not recovered reference families, accuracy
or a scientifically admitted result. Full raw evidence remains under
`benchmarks/results/qfo_sequence_graph_v1/all_hits/`.

## Independent Admission

Submitted21823 using the existing checked batch entrypoint after confirming
the terminal scheduler record and exact result hash. Frozen admission executor
`f1e21b09c28f270dc3ef2243bdcad86f212b58a0` has clean tracked scientific code.
Plan SHA-256 remains
`0e20d71936b92d7d2f02f3f35e7b7f380553803d9e482fde4d4b6c549b6e96d0`.
Admission was confirmed RUNNING with2CPUs/192GiB. It independently checks
source/checkpoint tuples, clustering evidence and output membership. Its
successful report, not this execution summary, is required before conversion.

Top100job21814 started after21813 and remains a separate diagnostic with
the same32CPU/384GiB allocation. No cap, scientific parameter or failure
policy changed based on the all-hit output.39focused admission/batch tests pass.

## Resource Scope

The retained [GNU-time record](qfo_sequence_all_hits_21813_time.txt) reports
32:50.27wall,1786.07user CPU seconds,183.86system CPU seconds,99%CPU and
40048236KiB maximum process RSS for the checked replay worker. The Slurm
37:34includes outer validation and preparation. Neither is end-to-end sequence
search timing. Allocating32CPUs does not imply32cores were used throughout.
Maximum process RSS is not simultaneous summed process-tree memory or a
cgroup peak. Shared-host load prevents controlled speed comparisons.

The previously malformed live sstat CPU value is not used. No missing
resource field is imputed. Pair conversion, six-endpoint scoring, independent
score admission and paired uncertainty all remain outstanding for this arm.
