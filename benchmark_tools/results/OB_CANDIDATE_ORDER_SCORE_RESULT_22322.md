# Candidate Hit-Order Effect

Job 22322 completed `0:0` in 1:14 using the separately prepared v2 plan.
The [submission](ob_candidate_order_score_submission_22322.json) and
[independent readback](ob_candidate_order_score_readback_22322.json) retain
the experiment. The [predecessor startup failure](OB_CANDIDATE_ORDER_FAILURE_22321.md)
remains recorded and was not resumed or relabeled as successful.

All five arms use the same 62,885 seed groups, gene indexing, species
memberships, frozen source, private runtime and satellite_v2 parameters.
The four factorial arms contain 18,235,373 identical unique directed nonself
hit keys; the separate full-fresh control adds 251,135 original self hits.
All five outputs cover 251,378 genes exactly once and contain 54,445 candidate
families. Independent reconstruction validates every candidate merge trace
and seed-family sidecar; all ten reported partition contrasts reproduce.

| Hit order | Normalized scores | Self hits | Genes changed versus historical candidates | Genes changed versus 0.11.0 fresh replay |
|---|---|---|---:|---:|
| Historical | Historical | No | 0 | 93 |
| Historical | Fresh | No | 0 | 93 |
| Fresh | Historical | No | 93 | 0 |
| Fresh | Fresh | No | 93 | 0 |
| Fresh | Fresh | Original self hits | 93 | 0 |

Thus changing hit order alone reproduces the entire residual candidate
partition discrepancy in this fixed-seed/runtime intervention. Changing only
the observed tiny score differences does not change either candidate partition.
Removing self hits also has no effect in this case. Both anchor outputs are
reproduced exactly as label-invariant partitions, not just equal group counts.

This is a tested candidate-stage order effect, not a claim that all candidate
inference is order-sensitive, a proof of a particular floating-point operation's
role, or final-F1 attribution. The unchanged implementation uses stable
cluster-pair aggregation and score-prioritized bounded attachment; numerical
aggregation/tie ordering is a plausible more detailed mechanism, but the
factorial does not isolate individual arithmetic operations.

Seventy focused tests pass across the external reader, driver alignment,
graph alignment and candidate auditors. The independent reader binds the
plan, source/runtime records, submitted command and terminal scheduler result,
rechecks all input/output identities, reconstructs candidate consistency and
compares partitions against both retained baselines. `env -i` intentionally
removes scheduler variables in the scientific child; job identity is linked
through the submission receipt and scheduler accounting, not an invented
child `SLURM_JOB_ID`. Shared-host duration is not controlled timing evidence.

No scientific defaults, dependency pins or historical benchmark scores were
changed. Next derive a small regression fixture from this mechanism and test
an explicitly versioned deterministic ordering policy without modifying the
frozen publication baseline. Any candidate-algorithm revision requires renewed
end-to-end validation; this result alone does not justify a new accuracy claim.

```bash
python -m benchmark_tools.audit_ob_candidate_order_scores --repo . \
  --directory benchmarks/work/ob_candidate_order_scores_v2_20260926 \
  --output /tmp/ob-candidate-order-readback.json
```
