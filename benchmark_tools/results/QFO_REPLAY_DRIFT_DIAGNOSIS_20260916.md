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

## Completed Exact-Entry Capture

Job21305 completed0:0 in9:42 using frozen observerf6da2bc and launcher49ab.
Before/after historical input, runtime and installed-package checks passed.
Independent verification confirmed saved graph fingerprints, source/input
records and complete partition comparison. Gene-name ordering is also
byte-identical to the checkpoint (SHA256
246d02da9635576f04e94a51dbec4093d33f3481bf7f761e03afb3a53d51a0d3).

All three initial RBNH arrays are byte-identical to diagnostic21295, but the
first clustering produced349,950 groups rather than349,898. There are4,427
diagnostic-only and4,479 capture-only groups; all976,504 genes are retained.
Capture partition SHA256:
def06c421941743617dd3540a14f8d9dc4f19eeedd34c1eb0c86c92160500488.
Its1,369,532 singleton edges differ from both diagnostic21295/historical
1,361,622 and replay21288's1,360,934.

This localizes the observed difference to the first clustering execution,
despite identical recorded input graph arrays and gene order, before profiles
or subsequent refinement. It does not yet identify the cause: the comparison
spans different processes/jobs, and complete native-library binary and
worker-environment identity across all historical runs is not established.
The observer adds instrumentation and an early stop. Do not treat this as
proof of a general algorithmic nondeterminism defect or a historical replay.

Evidence: `qfo_replay_initial_capture_20260916.json` and
`qfo_replay_initial_capture_comparison_20260916.json`.
Next run bounded fresh-worker repeats directly on the preserved identical
graph, recording worker inputs, command, native-library binaries and runtime
settings. Keep all partitions; do not select the repeat closest to historical
accuracy or restart full profile inference before this is understood.

## Three Preserved-Graph Repeats

Job21307 completed0:0 in32:26 using frozen executora69d505. Three sequential
fresh instrumented workers consumed the same preserved graph and gene order,
CPM0.1/seed4/include-isolates, with observed CPU affinity[8]. Each produced
349,898 groups and the identical partition SHA256
8c162782bba955078a8f7e4e6088df5b34c9f15b055803224e513fdea8b33cdd.
All match diagnostic21295 exactly and differ from capture21305.

Imported module files,71 loaded shared libraries, Python/package versions,
recorded environment, host/platform and affinity agree across these repeats.
Independently rechecked graph/source/native file hashes and every complete
partition comparison. Evidence: `qfo_saved_graph_repeats_20260916.json`.
The three worker durations were627.29,693.45 and555.00 seconds; these are
shared-node diagnostic costs, not a controlled performance comparison.

These three observations agreed under this instrumented single-CPU configuration.
They do not establish general determinism or identify the earlier discrepancy;
the later affinity panel below contains a same-affinity counterexample.
The matching diagnostic also requested one CPU, whereas the replay/capture
requested32; exact affinity and loaded-library records were not captured for
those earlier workers. CPU affinity and process/import context are candidates
for controlled tests, not established causes. Next vary one factor at a time
on the same preserved graph, starting with repeated one-versus32-CPU affinity
within a single allocation while keeping the worker code, inputs and settings
fixed. Do not change inference defaults or select the best partition by score.

## Completed Affinity Panel: Same-CPU Disagreement

Job21311 completed0:0 in43:09 from frozen executorcc3bffa. Four fresh workers
alternated one and32 available CPUs in the same allocation, with unchanged
graph/gene order, CPM0.1, seed4, include-isolates and one-thread OpenMP/BLAS
settings. No accuracy was evaluated and no output was selected or retried.

| Arm | Actual Available CPUs | Groups | Partition SHA256 Prefix |
| --- | ---: | ---: | --- |
| one_cpu_0 | 1 | 349,898 | 8c162782 |
| all_cpus_0 | 32 | 349,898 | 8c162782 |
| one_cpu_1 | 1 | 349,950 | def06c42 |
| all_cpus_1 | 32 | 349,898 | 8c162782 |

Both single-CPU workers used CPU8. All recorded software identities match
excluding the intentionally varied affinity, and affinity also matches within
each repeated arm. The second one-CPU partition has4,479 groups absent from
the first, while the first has4,427 absent from the second. Every partition
covers the same976,504 genes. Both32-CPU repeats are byte-identical. Neither
these matches nor the earlier three matches establish general determinism.

Independent admission re-read worker/execution records, checked248 unique
source/library/interpreter/graph/partition files, and recomputed all13 recorded
full-partition comparisons. Evidence:`qfo_affinity_verified_20260916.json`,
SHA256a927b00477f12e0a4c6548271942bc2f78d047bf0f7d0304d8aa0abd11954ec0.
The underlying job report SHA256 ise9311a5862b1db605b7d741003b2033de6f41f2c6db7e726316c181d00105de6.
Ten admission tests cover incomplete panels, changed identity, affinity,
graphs, arm order and exit status; full unit suite1,127 tests passed.

CPU availability alone cannot explain this panel's variation. Matching files
and selected runtime settings do not capture all native state, so the cause
is still unresolved. The installed Python binding passes the requested seed
to the optimizer; static inspection alone cannot prove its effective native
state or exclude graph-conversion differences. No library defect or historical
environment explanation is established by this experiment.

Next instrument the actual igraph edge order and weights at the optimizer
call boundary, alongside effective arguments and the resulting partition.
Compare bounded fresh-worker repeats without choosing by accuracy. This can
separate variability before the optimizer call from variability during native
optimization, but instrumentation itself must be documented. Do not substitute
a preferred partition, restart full QfO profiles, or declare the publication
baseline reproducible before stronger evidence is obtained.
