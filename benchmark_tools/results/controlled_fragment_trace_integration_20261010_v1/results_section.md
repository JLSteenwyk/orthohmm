### Retained Fragment Pipeline Observations

The table counts the fixed representatives (HS 12, PHY 18, OF 11,
OF checkpoint 12) in each arm, including errors, recoveries and retained
TP controls. These are not rates, prevalence estimates or new accuracy
endpoints. All 53 cases and empty bins remain available; no preferred
stage explanation was used in selection.

| Method | Recorded Decision | Baseline Representatives | Fragment Representatives |
| --- | --- | ---: | ---: |
| HS | group retention | 9 | 6 |
| HS | group separation | 3 | 6 |
| PHY | candidate separation | 1 | 4 |
| PHY | event rule retention | 6 | 7 |
| PHY | observed duplication exclusion | 6 | 4 |
| PHY | unambiguous bypass retention | 5 | 3 |
| OF | mcl group separation | 0 | 2 |
| OF | native pair retention | 8 | 7 |
| OF | within mcl native pair exclusion | 3 | 2 |
| OF checkpoint | mcl group separation | 3 | 3 |
| OF checkpoint | native pair retention | 9 | 9 |

Three illustrative records show distinct observable boundaries, not causal
or biological validation. HS Case0000 (unflagged FN) lacks direct
significant hits in either direction, yet its graph is connected and its
groups are separate. PHY Case0026 (two truncated endpoints, FN) has both
directed significant hits and connected graph support, but its observed
ancestor is called duplication and excludes the native pair. OF Case0030
(unflagged FN) is co-clustered at MCL but absent from native pairs; a
matching event adapter is unavailable, so its tree-level cause is unknown.

Search-hit absence does not distinguish prefilter rejection from scoring
or significance thresholds. Graph/group support is not orthology. A
recorded duplication establishes the pipeline's decision, not correctness
of an inferred tree or the underlying evolutionary history.

[Selection and empty bins](controlled_fragment_trace_selection_20261009_v1/selection.json),
[all retained observations](controlled_fragment_stage_trace_20261010_v1/report.json),
[stage table](controlled_fragment_stage_trace_20261010_v1/stages.tsv),
[independent readback](controlled_fragment_stage_readback_20261010_v2.json),
[actual execution and preserved refusal](controlled_fragment_trace_execution_20261010_v1.json).

