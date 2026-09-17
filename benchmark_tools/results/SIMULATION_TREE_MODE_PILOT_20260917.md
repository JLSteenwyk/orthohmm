# Supplied-Own-Tree Pilot

Dataset: baseline_20261101, the first baseline seed. This is an execution-mode
control, not a new accuracy comparison or evidence for general equivalence.
Both frozen tools reran from806original FASTA sequences with their own previously
inferred species tree supplied, four CPU threads and fresh output directories.

Job21332 completed0:0 in35s, including preflight/postflight checks. Executor
0e797bab52f91e494c1e3be4b728de3dbedd3fe8. The prior submission21331 failed its
commit-argument guard before inference; it is retained as a submission failure.

| Observation | OrthoHMM satellite_v2 | OrthoFinder3.1.5 full |
|---|---:|---:|
| Native output admission | Passed | Passed |
| Original native ortholog pairs | 2,881 | 2,884 |
| Supplied-own-tree native ortholog pairs | 2,881 | 2,884 |
| Added or removed native pairs | 0 | 0 |
| Rooted species-tree clades | Identical | Identical |
| Retained artifact bytes | All27identical | Differences retained below |

OrthoFinder retained201versus202 inventoried comparison artifacts. The omitted
file is Alignments_ids/SpeciesTreeAlignment.fa, associated with inference of the
species tree that this control supplies. Both MCL files differ in their command
comments' absolute paths. A separate reread using the audited read_checkpoint
parser and SequenceIDs mappings recovered98identical canonical groups covering
all806genes in both runs. The report retains these byte differences rather than
masking them. Other retained graphs, ID maps, gene-tree and alignment files match.

Snapshot simulation_mode_control_baseline_seed1_20260917.json SHA256
446a3c19cebf87fe562cc9e88cf183986e4e644b71116eaa5ba2839c9083f561.
Raw outputs: benchmarks/results/simulation_mode_control_baseline_seed1_v1.
No score against evolutionary truth was calculated. Shared-host wall time is not
a matched performance result. These comparisons do not establish equivalence
across the other69datasets or validate use of postprocessed output as a restart
checkpoint.

Independent admission completed in admit_simulation_mode_pilot.py: rechecked
the fixed scheduler job/executor, all native output inventories, original inputs
and baseline artifacts, reconstructed both commands independently, and recomputed
native-pair/topology/retained-artifact comparisons. The only permitted byte
exceptions are one MCL command comment per changed MCL file and the omitted
species-tree alignment; changed MCL matrix content is rejected.
Snapshot simulation_mode_pilot_verified_20260917.json SHA256
4a4630c8c036a4ae85045045ebf848ab9eefa730de9dcc8295570a50deeddd0c.
This authorizes expanding the controls, not assuming the remaining cells pass.
