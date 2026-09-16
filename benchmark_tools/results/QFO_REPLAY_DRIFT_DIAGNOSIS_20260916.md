# QfO Replay Drift: Initial Localization

The completed native replay21288 passed source/runtime/input integrity checks
but failed historical partition equivalence. It is not a successful historical
reproduction, and no new accuracy score has been computed.

## Earliest Recorded Difference

| Stage Count | Historical Production | Frozen Replay |
| --- | ---: | ---: |
| Significant hits | 88,729,858 | 88,729,858 |
| Initial RBNH edges | 24,148,515 | 24,148,515 |
| Initial singleton-assignment edges | 1,361,622 | 1,360,934 |
| Profiles built | 57,883 | 57,856 |
| Profile candidates | 19,519,955 | 19,519,973 |
| Significant profile hits | 1,966,961 | 1,963,669 |
| Strict profile edges | 98,085 | 98,795 |
| Final singleton-assignment edges | 1,250,874 | 1,236,855 |
| Final network edges | 25,168,957 | 25,159,378 |
| Final groups | 390,817 | 390,657 |

Source: historical `qfo_benchmark/results/orthohmm_high_sensitivity_isolated/metrics.json`
and `benchmarks/results/publication_qfo_replay_check_v1/replay.json`.
These are stage counts, not accuracy statistics. The first recorded difference
precedes profile expansion; profile scoring alone cannot explain that earlier
count difference. Equal initial edge counts do not prove identical graphs.

## Code Inspection

Historical694a77f and frozen7f3a9e4 have identical `externals.py` clustering
and `helpers.py` graph representation code. Both use Leiden CPM0.1/seed4,
include all input genes, and isolate graphs over five million edges in a
fresh worker. Both consume in-memory numeric edges; text-edge rounding is
not the path used for these calls.

The RBNH builder differs by an optional positive threshold factor whose
default is1.0. The next experiment tests the complete graph arrays rather
than assuming this change is neutral from inspection. The checkpoint
preserves hit array order and numeric values. Both profile branches use the
unrefined multipass clusters as seeds; the replay's separately written
refined checkpoint is not its profile-building input.

## Diagnostic Experiment

Run `diagnose_qfo_initial_graph.py` from a committed isolated executor:

- Verify the frozen launcher/runtime and historical checkpoint/metrics hashes.
- Load historical and frozen RBNH builders directly from their Git blobs,
  using verified identical graph helpers, and compare gene order plus exact
  source/target/weight array hashes.
- If graphs differ, stop before clustering and report that difference.
- Otherwise preserve the numeric graph and cluster it twice in separate
  workers with the same settings; record complete partitions and singleton
  edge arrays/counts, and compare the partitions exactly.
- Do not rerun sequence searches, profiles, reference scoring, or tune any
  inference settings based on benchmark outcomes.

Historical intermediate partitions, graph hashes and complete historical
native/package inventories are unavailable in the retained run records.
The experiment can test current graph-source equivalence and repeatability;
it cannot by itself certify or reconstruct the historical runtime. No cause
is established yet, and no automatic restart of the full replay is planned.

## Completed Initial-Graph Experiment

Job21295 completed0:0 in24:18. Historical and frozen builders, both evaluated
in the current environment on the same checkpoint, produced byte-identical
source/target/weight arrays for24,148,515 RBNH edges. This rules out the
threshold-factor source change as an explanation for this graph. It does
not compare against an archived historical graph, which is unavailable.

Two separate initial-clustering workers produced byte-identical349,898-group
partitions (SHA2568c162782bba955078a8f7e4e6088df5b34c9f15b055803224e513fdea8b33cdd).
Both yielded1,361,622 singleton-assignment edges, matching the historical
recorded count rather than the failed-equivalence replay's1,360,934.
All singleton source/target/weight fingerprints also agree between repeats.
Recorded environment: NumPy2.2.6, igraph1.0.0, leidenalg0.11.0.

Independently rechecked saved source/checkpoint/graph-file hashes, graph-array
fingerprints, complete input-gene coverage and exact partition comparison.
Evidence: `qfo_initial_graph_diagnostic_20260916.json`.

The discrepancy is not resolved. In particular, the failed-equivalence
replay did not retain its initial RBNH arrays or initial partition. Next
instrument that exact replay entry point to capture these before proceeding
to profiles, then compare them with the preserved deterministic diagnostic.
Do not infer general Leiden nondeterminism, an HMM-profile defect, or full
historical reproduction from the current evidence. No accuracy was scored.
