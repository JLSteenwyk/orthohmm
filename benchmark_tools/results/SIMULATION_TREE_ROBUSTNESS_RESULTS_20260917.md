# Simulation Tree Robustness Results

The prespecified complete panel is now scored. These are exploratory,
development-exposed, simplified simulations, not independent biological
validation or matched runtime measurements. No settings or seeds were selected
from accuracy outcomes.

## Admission and Provenance

Independent native auditor21435 completed0:0 in4:29. Its420-method inventory
admits210/210 OrthoHMM supplied-tree runs and195/210 OrthoFinder runs. The15
OrthoFinder failures are divergent seeds20261101,20261102,20261107,20261108,
20261109 under each of the three trees: nonfinite native graph weights. These
same five original inferred-tree baselines failed. Supplying a tree does not
repair this upstream problem. All requested topologies were retained in the
405 admitted outputs.

The original OrthoHMM inferred-tree arm has three failures, while all supplied
arms complete. This is a completion observation, not a measured accuracy gain
for those missing inferred baselines. The analysis retains560 total arm rows:
537 scored and23 failed, with no failed-run accuracy imputation.

Raw admission remains outside Git at
`benchmarks/results/simulation_tree_panel_admission_v1/results.json`:
SHA256 `17a1be71823b9fb7fb13082d05001d7667bb89c3366da2e3138a124e90907820`.
It contains the complete native inventories; archival packaging remains due.

The [cross-arm artifact audit](simulation_tree_artifacts_verified_20260917.json)
has400 retained-upstream-equivalent,2 different and18 unavailable contrasts.
SHA256 `047dd998752dd30ab62e042e6b3e5f3442ec4f4dffcff030984f3e01e4fff65c`.
Both differences involve OrthoFinder divergent seed20261106: generating versus
inferred, and NNI1 versus generating. OG0000004 alignment and raw gene-tree
files differ beyond permitted MCL command-path comments and species-alignment
availability. These cases were not excluded. They prevent a strict tree-only
causal attribution for those comparisons; retained-file equivalence elsewhere
does not prove identity of all internal computation.

Scoring used frozen worktree `publication_simulation_tree_scoring_v1`, commit
04055a5. It rechecked inputs, truth and native prediction hashes and reproduced
every successful original inferred score exactly. Raw complete scores:
`benchmarks/results/simulation_tree_scores_v1/results.json`, SHA256
`58df72a96ed0686380b2e9a5c9bbffc7087e4985f9b8a809393f137586048b8f`.
The [complete summary](simulation_tree_robustness_summary_20260917.json) has
SHA256 `1f965d18f572864922e1cdaffb9e4f0443d4ba91ab8b7418dbc99482a162674b`.

## Observed F1 Effects

Values are mean paired-seed F1 differences in percentage points, not pooled
pair-count scores. NNI contrasts use the generating-tree arm as reference.
The supplied generating tree is an oracle diagnostic, not an achievable
end-to-end baseline or a guaranteed upper accuracy bound.

| Condition | Method | Generating - inferred | NNI1 - generating | NNI2 - generating |
| --- | --- | ---: | ---: | ---: |
| Baseline | OrthoHMM | 0.035 | -0.282 | -0.410 |
| Baseline | OrthoFinder | -0.005 | -0.198 | -0.137 |
| Divergent | OrthoHMM | 0.000 | -0.134 | -0.136 |
| Divergent | OrthoFinder | 0.000 | 0.000 | 0.000 |
| Turnover | OrthoHMM | 0.148 | -1.207 | -2.341 |
| Turnover | OrthoFinder | 0.019 | -1.524 | -1.604 |
| Divergent turnover | OrthoHMM | 0.127 | -0.713 | -1.333 |
| Divergent turnover | OrthoFinder | 0.017 | -0.225 | -0.494 |
| Missing20 | OrthoHMM | 0.008 | -0.260 | -0.456 |
| Missing20 | OrthoFinder | 0.160 | -0.149 | -0.211 |
| Uneven taxa | OrthoHMM | 0.050 | -0.296 | -0.786 |
| Uneven taxa | OrthoFinder | 0.086 | -0.174 | -0.690 |
| Taxon-count control | OrthoHMM | 0.217 | -0.267 | -0.727 |
| Taxon-count control | OrthoFinder | 0.085 | -0.352 | -0.870 |

OrthoHMM generating-versus-inferred uses9 divergent and8 divergent-turnover
seeds; its other comparisons use10. OrthoFinder uses5 divergent seeds and10
elsewhere. Complete-case results may be biased by failures.

The fixed analysis uses20,000 paired-seed bootstrap replicates, seed20260918,
and Bonferroni adjustment across all126 planned F1/precision/recall endpoints.
None of the generating-versus-inferred adjusted intervals excludes zero.
NNI2-minus-generating F1 and recall intervals are below zero for OrthoHMM
turnover, divergent turnover, uneven taxa and taxon-count control; and for
OrthoFinder uneven taxa and taxon-count control. These are12 adjusted endpoints
in total. No precision interval excludes zero. All intervals, seed effects,
exclusions and undefined-ratio flags are retained in the machine-readable
summary; a zero-width bootstrap interval is not proof of biological equivalence.

Ten planned seeds provide limited tail resolution, especially for126-way
adjustment. Percentile intervals have approximate, not guaranteed simultaneous
coverage. The perturbations are deterministic rooted NNI choices with branch
lengths carried along subtrees, not sampled empirical tree uncertainty. No
claim of arbitrary-tree robustness, comparative superiority, real-domain or
fragment behavior follows.

## Figure and Manuscript

The [complete nine-panel figure](figures_simulation_tree_robustness_20260917/simulation_tree_robustness.png)
shows all126 endpoints, with nominal and adjusted intervals, method-specific
symbols and shared scales down contrast columns. PDF/SVG versions and a
source/output hash manifest accompany the PNG. The plotter recomputes the
fixed summary from all560 records and rejects altered contrasts or bootstrap
settings before rendering. It distinguishes unavailable effects and single-pair
descriptive effects from estimated intervals. The caption and manuscript retain
exact paired-seed counts, upstream differences and failure-conditioning caveats.
Methods and Results are integrated into the publication manuscript draft.
