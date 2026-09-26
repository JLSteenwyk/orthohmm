# Matched-Recall Initial-Search Graph Control

## Frozen Scope

Follow the completed search-only matching gate, without changing its selected
DIAMOND cutoff (reported E <= 1e-40) or fixed HMM search (native E < 1e-4).
Use only the **35 reporting datasets**, seeds 20261106-20261110 in all seven
conditions. The 35 calibration datasets are excluded from this downstream
comparison. Freeze this protocol before generating either new graph partition.
All datasets remain development-exposed, including earlier full-method results;
this is a controlled diagnostic, not independent confirmation.

Pin the admitted search report
`search_sensitivity_results_20260926.json`, SHA-256
`d58b1603e563eda27d90cc113950fb28091c94d57233a60ec27b32e297ad729c`,
and its execution receipts/raw hits. Reuse these searches, never silently rerun
or rescore sequences. No cutoff, candidate cap, clustering resolution, seed,
refinement parameter or endpoint may be selected using downstream scores.

## Two Arms and Shared Inference

Both arms use the installed frozen OrthoHMM scientific source and runtime from
`publication_frozen_overlay_install_20260926.json`. Call the installed production
graph APIs in isolated Python, not an unpinned development checkout. Use
lexically sorted gene IDs and species filenames for deterministic indices.
Retain within-species and self hits as production graph inputs; only evaluation
excludes within-species/self pairs. Include genes with no retained hits.

- **HMM:** reuse all retained native hits and their already length-normalized
  scores unchanged. Do not divide by sequence length again.
- **DIAMOND:** retain hits with reported E <= 1e-40, without another target cap.
  Use raw score divided by sqrt(query length * target length), matching the
  existing sequence-search control adapter. Do not substitute bit scores.

Reject duplicate directed pairs, foreign IDs, inconsistent lengths, nonfinite
scores or malformed records. Preserve exact numeric graph inputs/checksums.
Matching recall does not make the two score distributions, hit identities or
rankings equivalent; this experiment changes the search evidence as a whole.

For both arms, execute the existing profile-off graph path:

1. `build_rbnh_edges` with all indexed hits and species assignments.
2. Leiden CPM resolution 0.1, seed 4, including isolates.
3. `build_singleton_assignment_edges`, combine with the original edges, and
   repeat Leiden with the same settings.
4. Production `refine_cluster_indices` with the same graph/hit inputs and
   production defaults. Respect the existing broad-panel copy-only input rule,
   although these datasets have fewer than its ten-species threshold.

Profile expansion, candidate-family expansion and phylogenetic reconciliation
are **off in both arms**, as in the previous initial-search controls. The final
refined partition is the sole scored output; intermediate partitions/graphs
are diagnostic artifacts, not alternative endpoints. Preserve all genes once,
including singleton groups. No truth file enters inference.

## Evaluation and Uncertainty

Convert each final partition to canonical unordered **cross-species pairs**.
Score against the retained simulator's speciation-derived `ortholog_pairs`,
not ancestral-family membership used for search recall. Reuse the validated
simulation pair scorer, independently checking input coverage, pair identities
and TP/FP/FN arithmetic. These cluster-derived pairs are not native reconciled
ortholog predictions. Report gene coverage, nonsingleton coverage and group
counts alongside pair precision, recall and F1.

Primary contrast: HMM minus DIAMOND mean seed-level pair F1, averaged equally
across the seven conditions. Report all seven condition contrasts too.
Precision/recall and coverage are explanatory, not alternative winner-selection
criteria. Keep per-seed effects and wins/ties/losses.

Use paired bootstrap of the five reporting **seed blocks**, carrying all seven
conditions together per sampled seed. Recompute condition means and the overall
equal-condition mean per replicate. Fix PCG64 seed 20260927 and 20,000
replicates. Report marginal 95% intervals and Bonferroni-adjusted percentile
intervals for the eight F1 contrasts (overall plus seven conditions).
Precision/recall intervals, if displayed, are exploratory marginal intervals.
Five blocks give limited resolution and approximate coverage; do not treat
gene pairs, families or shared-history conditions as independent replicates or
claim population-wide significance from these intervals.

Require both arms on all 35 datasets for complete-panel contrasts. Preserve
failures as missing with their causes, never as zero accuracy. Do not retry
native failures automatically, drop a difficult condition or switch to an
intermediate partition. A failure prevents the complete-panel claim.

## Execution and Claims

Use immutable executor code, one attempt per cell, private fresh output paths,
before/after source and input hashes, exact command/environment receipts,
per-stage logs/resources and terminal scheduler verification. Up to two
4-CPU/8-GiB jobs may execute concurrently on the shared host. Timings measure
incremental conversion/graph work and are descriptive, not dedicated matched
efficiency evidence. Preserve main OrthoBench job 22179 and all prior results.

The result may support only a statement about initial-search evidence under
these settings and matched simulation recall. It does not isolate profile
expansion, assess the full phylogenetic pipeline, establish real-data matching,
or justify transferring E <= 1e-40 to QfO/OrthoBench. No scientific defaults
change. Retain negative and neutral outcomes and all broader publication gaps.
