# Downstream Dependency Replay

Job 22320 completed with scheduler exit `0:0` after 20:40. Both predeclared
arms succeeded once, with no retries. The
[independent readback](ob_dependency_replay_readback_22320.json) checks the
pinned plan, recorded commands/environment, scientific source and native
extension identities, four graph snapshots, retained partitions, and candidate
seed/merge consistency. All partitions cover the same 251,378 genes exactly
once. No accuracy evaluation or phylogenetic inference was run.

Both arms use the same fresh installed search checkpoint and unchanged
scientific code, CPM 0.1, seed 4, one BLOSUM62 profile iteration, and satellite_v2
candidate expansion. The contrast changes the complete retained leidenalg
distribution, including bundled native libraries, not just a version label.

| Stage | Distribution 0.12.0 groups | Distribution 0.11.0 groups | Genes in changed groups |
|---|---:|---:|---:|
| Initial clustering | 58,827 | 58,772 | 23,352 |
| Multipass clustering | 51,130 | 51,181 | 26,935 |
| Multipass refined | 63,227 | 63,245 | 19,669 |
| Profile-base clustering | 57,742 | 57,696 | 29,745 |
| Profile-final clustering | 50,857 | 50,894 | 33,061 |
| Profile-final refined | 62,928 | 62,885 | 25,335 |
| Candidate families | 54,495 | 54,445 | 30,788 |

The initial graph arrays are exactly equal: 1,803,122 edges in both arms.
Subsequent graph endpoints differ, consistent with propagation through
group-dependent processing. The candidate traces reconstruct 8,433 and 8,440
merges respectively and agree with their seed-family sidecars. This checks
recorded merge consistency, not independently recomputed search support.

The 0.12.0 candidate partition exactly matches the fresh installed full run.
The 0.11.0 result shares 54,433 candidate groups with the historical partition;
12 groups on each side differ, involving 93 genes. Equal total group counts
therefore do not establish exact historical reproduction. Both historical
comparisons use partition records from the pinned prior search-comparison audit.

This establishes a dependency effect extending through candidate formation.
It does not explain the residual 93-gene discrepancy, isolate an upstream
library component, or attribute the final 0.2845-point F1 decrease. The latter
requires downstream phylogenetic evidence. No default or dependency pin was
changed, and no historical score was replaced. One run per arm does not
establish general determinism. Arm durations of 656.883 and 580.825 seconds
are shared-host diagnostic observations, not comparative timing results.

One hundred focused tests pass across the new reader, candidate auditor,
graph checks, replay driver, search comparison and installed-run reader.
The first local readback is retained; the second adds candidate comparisons
against historical/fresh baselines without rerunning inference. The committed
report is the latter readback.

```bash
python -m benchmark_tools.audit_ob_dependency_replay --repo . \
  --directory benchmarks/work/ob_dependency_replay_v2_20260926 \
  --output /tmp/ob-dependency-readback.json
```

Next locate the residual difference in retained historical intermediate
partitions before deciding whether another bounded native run is justified.
Publication-wide provenance, independent validation limitations, uncertainty,
controlled timing and release work remain open. DGX stays deferred.
