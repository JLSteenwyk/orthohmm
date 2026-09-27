# Matched-Graph Execution Binding Audit

The original readback checked hashes, numeric mappings, edge derivations,
partitions and commands, but did not enforce exact recorded environment or
bind every receipt inventory to the local files being evaluated. Membership
checks alone also allowed empty or duplicated native module lists.

The strengthened `audit_matched_graph.py` now requires:

- The exact runner environment, including single-thread library settings,
  fixed hash seed, private HOME and no additional variables.
- Private numeric input, native receipt and both stage logs at their expected
  locations with matching checksums.
- Native completion status, installed interpreter/prefix and the three
  required installed production modules, without missing or duplicated entries.
- All five expected graph/partition outputs and the expected checkpoint
  manifest, bound to the artifacts independently read and scored.

All 70 retained cells passed the
[additional readback](matched_graph_readback_v2_20260926.json). Its `cells`
array is exactly equal to the original readback: no inference outputs, evidence
records, scoring inputs, thresholds or scientific defaults changed. The
original report remains preserved and remains the report pinned by the scoring
artifact. This is an additive audit, not a replacement experiment.

The focused audit, scorer and runner suite passed 45 tests, including rejection
of changed/extra environment variables, redirected artifacts, incomplete logs,
incorrect executable/prefix and missing/duplicated/out-of-prefix modules.
The broader matched-graph, statistics, figure and resource suite passed 90 tests.

Reproduce from the retained local artifacts with a new output filename:

```bash
python -m benchmark_tools.audit_matched_graph \
  --submission benchmark_tools/results/matched_graph_submission_20260926.json \
  --output /tmp/matched_graph_readback.json
```

The environment is verified against the historical runner receipt, not an
independent in-process observation. This does not rerun Leiden, independently
derive refinement decisions, establish controlled timing, or supply independent
biological validation. All original scientific claim boundaries remain.
