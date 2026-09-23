# Publication Claim-To-Evidence Checklist

Status updated 19 September 2026. This is a completion audit, not a replacement
for the original publication goal. A linked plan or passing unit test is not
evidence that an experiment completed or a biological hypothesis is true.

## Claim Boundaries

| Proposed statement | Evidence | Assessment |
| --- | --- | --- |
| Satellite_v2 has a higher observed OrthoBench aggregate F1 than full OrthoFinder | [Paired analysis](ORTHOBENCH_UNCERTAINTY_20260916.md) | Descriptively supported; interval includes zero; development-exposed |
| Satellite_v2 trades higher precision for lower recall on OrthoBench | [Paired analysis](ORTHOBENCH_UNCERTAINTY_20260916.md) | Supported within this benchmark; not selection-adjusted generalization |
| OrthoHMM outperforms full OrthoFinder overall | [Eight-method comparison](PUBLICATION_COMPARISON_ORTHOMCL_COMPLETE_20260916.md) | Not supported; endpoints and benchmark rankings differ |
| HMM expansion contributes in historical OrthoBench processing | [Historical component audit](HISTORICAL_PROFILE_ABLATION_20260916.md) | Descriptive +0.595610 F1 points; current factorial intervals include zero |
| Initial HMM search improves F1 over sequence-search replacement | [OrthoBench controls](OB_SEQUENCE_SEARCH_RESULTS_20260916.md), [corrected QfO controls](QFO_SEQUENCE_UNCERTAINTY_RESULT_20260918.md) | Not established: adjusted F1 intervals include zero on both benchmarks. Corrected SwissTrees precision differences favor HMM; recall intervals include zero. Development-exposed evidence with unmatched search sensitivity/calibration/cost |
| Broad candidates improve reconciliation | [OrthoBench factorial](ORTHOBENCH_FACTORIAL_INTERPRETATION_20260916.md), [original-release QfO factorial](QFO_FACTORIAL_SWISS_RESULTS_20260918.md), [corrected factorial](QFO_CORRECTED_FACTORIAL_COMPLETE_20260919.md) | Candidate-expansion simple-effect F1 intervals include zero. Corrected SwissTrees C-by-R F1 interactions are positive with adjusted intervals excluding zero at both profile settings; original-QfO interaction intervals include zero. Conditional end-to-end interaction, not an isolated mechanism or general superiority |
| The frozen method was evaluated on novel taxa | [YGOB evaluation](YGOB_FROZEN_INTERPRETATION_20260916.md) | Complete bounded transfer evaluation; satellite F1 difference interval includes zero, precision higher and recall lower; not family-disjoint or superiority evidence |
| Validation is family-disjoint | [Homology screen](ygob_homology_screen_20260916.json) | Not established; substantial detected overlap |
| Original-release OrthoMCL final-group QfO scoring is complete | [Verified result snapshot](publication_comparison_orthomcl_complete_20260916.json) | Supported for the original inputs only; corrected-release BLAST 21713 was verified running on September 20 and no corrected OrthoMCL score is available |
| OrthoMCL BLAST failures have negligible impact | [Failure-impact audit](ORTHOMCL_FAILURE_IMPACT_20260916.md) | Not established; direct exposure is measured, indirect and counterfactual effects are not |
| Three Kingdoms demonstrates proteome-wide accuracy | [Supplementary score record](three_kingdoms_parity_20260907.json) | Unsupported; restricted BUSCO-reference universe |
| Every historical Three Kingdoms method used identical input bytes | [Method-input audit](THREE_KINGDOMS_METHOD_INPUTS_20260918.md), [completed matched rerun](THREE_KINGDOMS_SONIC_MATCHED_RESULT_20260918.md) | Not established historically: SonicParanoid native snapshot matches raw Danio rather than the staged stop-marker-stripped version; older high-sensitivity record lacks per-file hashes. Matched inference21795 and assessment21796 completed successfully, admitting contemporary Sonic F1=0.9912758996728462; this does not repair historical provenance or isolate the mismatch's causal effect |
| Historical Three Kingdoms scores reproduce from normalized groups | [Independent arithmetic audit](THREE_KINGDOMS_PAIR_COUNT_AUDIT_20260918.md) | Supported for all eight retained methods; native Sonic group conversion also verified. This does not establish matched historical inputs or proteome-wide accuracy |
| OrthoHMM is faster or more memory efficient under matched conditions | [DGX disposition](DGX_POSTRUN_ADMISSION_20260918.md), [descriptive observations](dgx_descriptive_resources_20260918.json) | Not established; all27 native-valid runs retained as descriptive evidence, with host-isolation uncertainty and asymmetric process-sampling gaps. No controlled comparison admitted |
| Outer PATH records prove OrthoFinder's historical companion-tool versions | [Child-PATH audit](DGX_SCALING_MIGRATION_20260917.md) | Unsupported: installed OrthoFinder rewrites its subprocess PATH; current reconstruction resolves bundled DIAMOND2.0.13/FastTree2.1.11/MCL14-137 instead of outer versions. Historical exec-path evidence still requires audit |
| ARM and x86 scoring are generally equivalent | [Portability diagnostic](native_scoring_portability_20260917.json), [development fix](NATIVE_BANDING_FIX_20260917.md), [one pipeline fixture](dgx_orthohmm_pipeline_smoke_20260917.json) | Not established: frozen narrow-band discrepancies remain; tested correction is in development source only, not the baseline; default64 matches tested synthetic fixtures and both OH modes match one simulation fixture only |
| A nearby parameter choice improves frozen-method OrthoBench F1 | [Six-variant panel](OB_PARAMETER_NEIGHBORHOOD_RESULTS_20260916.md) | Not established: all six adjusted F1 intervals include zero; no default promotion |
| The QfO native graph is reproducibly constructed | [Checked initial-graph repeats](QFO_CHECKED_REPEAT_RESULTS_20260917.md) | Three checked Python-pair runs preserve the full graph and yield identical partitions; not general determinism or complete historical replay equivalence |
| Profile-branch processing improves recovered-stage SwissTrees accuracy | [Paired SwissTrees intervals](QFO_SWISS_INTERVALS_20260917.md) | Not established: both profile contrasts have negative observed F1 differences and adjusted intervals including zero; effects occur in CASP and GH14 only. No superiority or equivalence claim |
| Recovered-stage QfO uncertainty is fully characterized | [Count audit](qfo_swiss_counts_20260917.json), [paired intervals](qfo_swiss_intervals_20260917.json) | Only SwissTrees completed: 18-family paired resampling, all12 adjusted intervals include zero. Other challenges and secondary mean require separate methods |
| OrthoHMM phylogeny exceeds full OrthoFinder on SwissTrees F1 | [Original-release intervals](QFO_SWISS_COMPARATOR_INTERVALS_20260917.md), [corrected-release intervals](CORRECTED_SWISS_COMPARISON_RESULT_21987.md) | Unsupported: the original-release adjusted F1 interval is negative; the corrected difference is -0.014900 with adjusted interval [-0.087078, 0.072090]. Neither superiority nor equivalence is established. Conditional evidence from 18 development-exposed families; no release interaction test |
| Domain architecture explains the OrthoHMM configuration effect | [Stratified results](SWISS_DOMAIN_STRATA_RESULTS_20260917.md), [all27 endpoints plotted](SWISS_DOMAIN_STRATA_FIGURE_20260917.md) | Not established: all nine adjusted interaction intervals include zero; annotation-defined association does not establish causality |
| Satellite constraints explain the five focal WGD homolog-coverage losses | [Prespecified case trace](BIOLOGICAL_WGD_CASE_TRACE_20260917.md) | Contradicted for these five homologs: all were in anchor candidates and separated during root-lineage reconstruction before constraints; topology correctness and upstream effects remain unresolved |
| The original-release QfO factorial has paired uncertainty estimates | [Validated results](QFO_FACTORIAL_SWISS_RESULTS_20260918.md), [counts](qfo_factorial_swiss_counts_20260918.json), [42 endpoints](qfo_factorial_swiss_bootstrap_20260918.json) | Complete for SwissTrees only: all adjusted F1 intervals include zero; R increases precision and lowers recall. No adjusted C-by-R interaction established for this release. Corrected-release results are complete and reported separately |
| Native GO/EC/FAS error bars establish paired method differences | [GO/EC audit](QFO_GO_EC_ARITHMETIC_AUDIT_20260917.md), [FAS sample audit](QFO_FAS_SAMPLE_AUDIT_20260917.md) | Unsupported: GO/EC use Student-t95% half-widths, FAS uses sample SEM, and none supplies dependency-aware paired method intervals |
| Supplying the generating tree improves simulation F1 | [Complete tree panel](SIMULATION_TREE_ROBUSTNESS_RESULTS_20260917.md) | Not established: all adjusted generating-versus-inferred intervals include zero; supplied-tree completion rescues three OrthoHMM baselines but does not supply their missing inferred accuracy |
| OrthoHMM is insensitive to species-tree error | [Complete tree panel](SIMULATION_TREE_ROBUSTNESS_RESULTS_20260917.md) | Unsupported: NNI2 F1 and recall deficits have adjusted intervals below zero in four conditions; bounded exploratory result, not arbitrary-tree robustness |
| The package is publication-ready | All sections below | Not achieved |
| Original QfO inputs match the corrected2020benchmark release | [Corrected archive comparison](QFO_CORRECTED_ARCHIVE_ACQUIRED_20260918.md) | Contradicted for the Xenopus proteome; preserve original results as release-limited |
| Corrected QfO inputs cover the retained reference identities and sequence content | [Native sequence and staging audits](QFO_CORRECTED_INPUTS_STAGED_20260918.md) | Supported: all 984,137 identities, 983,959 exact sequences and 178 representation-only differences; no unexplained differences. This is input compatibility, not biological annotation validation |
| Corrected-input cached replay reproduces the native high-sensitivity partition | [Independent replay admission](QFO_CORRECTED_REPLAY_ADMITTED_20260918.md) | Supported for this frozen run: 391,908 identical final groups across 984,137 genes, with 419 provenance records checked. Not historical-input equivalence, general determinism, accuracy or comparative timing |
| Corrected QfO accuracy or rankings are established | [Admitted partial table](qfo_corrected_comparison_20260923_v6/scores.md), [seven-method comparator intervals](qfo_fastoma_swiss_uncertainty_22098.json), [complete factorial](QFO_CORRECTED_FACTORIAL_COMPLETE_20260919.md) | Seven point-estimate rows and seven conditional SwissTrees contrasts are admitted; OrthoMCL remains unavailable. Full OrthoFinder has higher point estimates than phylogenetic OrthoHMM on all three F1 endpoints; OrthoHMM has higher GO/EC similarity and FAS. Corrected phylogenetic-minus-sensitive SwissTrees F1 is +0.148015 with adjusted interval [0.049791, 0.266015], but phylogenetic-minus-OrthoFinder intervals include zero. FastOMA uses a supplied species tree; its F1 difference versus full OrthoFinder is -0.068194 with adjusted interval [-0.158907, 0.026407]. All eight factorial cells are separately complete. No general superiority, equivalence or complete ranking; historical results cannot be relabeled |
| Original TreeFam family-level uncertainty can be recovered from pooled pairs | [Source retrieval investigation](TREEFAM_SOURCE_RETRIEVAL_20260918.md) | Unsupported: original trees and mapping remain missing; downloaded pooled reference is not an independent-family inventory |
| Merging overlapping VGNC labels makes ordinary family resampling valid | [Dependency audit](VGNC_DEPENDENCY_STRUCTURE_20260918.md) | Not established: 16,863 labels form 16,844 reference blocks, but almost all scored false positives cross blocks. Outcome-defined prediction components are not independent reference units |
| Missing initial OrthoBench edges explain final grouping errors | [Initial-edge trace](OB_INITIAL_EDGE_TRACE_20260918.md) | Descriptive localization only: 505 hit-supported separated memberships had initial edges and 1,466 did not. Other graph paths and later stages prevent a causal conclusion |
| Observed search rejections identify the cause of final OrthoBench F1 losses | [Raw search-decision recount](OB_SEARCH_DECISION_RESULT_20260918.md), [final-grouping join](OB_SEARCH_GROUPING_JOIN_20260918.md) | Rejection stages are localized on the watched set, not causally attributed: 48,368 directed pairs were prefilter-excluded and 1,619 scored non-significant, with zero historical accepted-presence disagreements. Yet 8,386 pairs excluded in both directions share final groups. No counterfactual scoring or authenticated historical raw-score equality |
| OrthoBench factorial statistics reproduce outside the checkout | [Relocated reproduction](ORTHOBENCH_FACTORIAL_REPRODUCTION_20260918.md) | Exact agreement for the eight-cell, 70-family, 20,000-draw statistical analysis; not native inference/scoring or cross-platform reproduction |
| Corrected SwissTrees sequence-control statistics reproduce outside the checkout | [Isolated and outside-repository reproduction](QFO_SEQUENCE_REPRODUCTION_20260918.md) | Exact agreement for all numerical fields of the three-arm, 18-family, 100,000-draw analysis from pinned audited counts; not a new source admission, native inference/scoring or cross-platform validation |
| Corrected SwissTrees factorial statistics reproduce outside the checkout | [Isolated reproduction](QFO_CORRECTED_FACTORIAL_REPRODUCTION_20260919.md) | Exact agreement for all numerical fields of the eight-cell, 18-family, 100,000-draw analysis and all 42 endpoints. Reproduces retained counts, not native inference/scoring, raw admission or cross-platform validation |
| New counter controls establish controlled comparative timing | [Earlier interval integration](DGX_INTERVAL_NATIVE_RESULT_20260918.md), [quiet hierarchy result](DGX_HIERARCHY_QUIET_RESULT_20260918.md), [native frontier result](DGX_FRONTIER_NATIVE_RESULT_20260918.md), [completed overhead audit](DGX_FRONTIER_OVERHEAD_RESULT_21838.md) | Not established: all earlier flags remain retained. Array21838 is terminal with six validated tasks, three wrapper failures and nine missing detailed scheduler records. Only one of nine paired comparisons is available, so the complete-panel overhead budget is unestablished. Transient service cgroups explain two identity-check failures, not their runtime effects. Non-CPU isolation and scientific inclusion rules remain unresolved. No historical timing is upgraded |

## Completion Requirements

The [completed lineage-native diagnostic](LINEAGE_NATIVE_RESULT_21995.md)
adds three validated, output-equivalent executions without the earlier
sibling-inventory failure. It retains 1/7/0 narrow flags in high-sensitivity,
satellite_v2 and full-OrthoFinder runs, respectively. This supports bounded
collector operation and output reproduction, not absence of interference,
overhead-budget passage or controlled comparative timing. The complete
[new paired overhead experiment](LINEAGE_OVERHEAD_RESULT_21999.md) now has
18 validated tasks and nine equivalent-work/duration pairs. All numerical
budgets pass, but 143 original and 23 narrow interval flags remain. Boundary
interval coverage is unavailable. No historical timing is upgraded.

The [complete interval description](LINEAGE_OVERHEAD_FLAGS_21999.md)
retains all 5,929 periodic intervals and all 23 narrow flags. Positive
root-minus-system.slice differences in flagged intervals do not identify
outside processes or establish causal interference. The
[read-crossing control](LINEAGE_READ_CROSSING_RESULT_22018.md) passed all
three prespecified finite-service trials with independent raw replay.
This establishes the tested lifecycle/read-window behavior only; it does
not explain native flags or establish root/user-slice specificity, general
churn robustness, non-CPU isolation or scientific timing eligibility.

The [root-context panel 22019](ROOT_CONTEXT_RESULT_22019.md) validates only
its first three native-only trials. Its first user-service condition failed
after the login user manager shut down; eight remaining conditions were
unrun. No intended-load sensitivity or timing-validity claim follows.

The separately retained [held-session follow-up 22020](ROOT_CONTEXT_RESULT_22020.md)
validates all 12 controls and the three prespecified aggregate user-load
responses. Common-work flags are 57/57 in user-contended and 0/171 in
idle/native-only conditions. A restarting unrelated daemon and the waiting
session remain documented parts of the environment. This supports only the
tested monitor response, not prior-flag attribution, zero observer overhead,
general isolation or scientific comparative timing admission.

The [native root-context integration 22021](ROOT_CONTEXT_NATIVE_RESULT_22021.md)
validates all three native commands, raw replay and canonical equality to
prior lineage outputs. It retains 52 original/11 narrow flags in 1,976
intervals and documents 383 scheduled restarts of an unrelated user service.
The claim is native integration only: no overhead bound, causal attribution,
background-free execution or scientific timing admission is established.

The subsequent [paired overhead panel 22022](ROOT_CONTEXT_OVERHEAD_RESULT_22022.md)
validates all 18 tasks and nine pairs against the native baseline. Every pair
and complete method median passes the prospective +10%/+5% engineering
budgets. All 344 original and 60 narrow flags remain. This supports the bounded
incremental collector comparison, not causal overhead, background-free
execution or a native-method speed ranking. Scientific timings remain unadmitted.

The [non-CPU environmental assessment](DGX_NON_CPU_ASSESSMENT_20260920.md)
retains positive native-step I/O stalls in all 18 tasks, memory stalls in 12,
and zero recorded memory-limit/OOM events. These are observed native-step
conditions, not identified outside interference or evidence of thermal/GPU/
memory-bandwidth isolation. No pressure-based exclusion rule is introduced.

The completed [full-node controls](FULL_NODE_CONTROL_RESULT_21918_20260919.md)
validate all nine injected workloads: zero narrow flags in 126 native-only
intervals and detection in 63/63 contended intervals. This supports bounded
monitor response, not controlled comparative timing. [Pressure and residual
descriptions](FULL_NODE_CONTROL_DESCRIPTION_20260919.md) show that process
creation and known contention both increase native pressure; pressure alone
does not identify outside work. [Longer-window diagnostics](CPU_WINDOW_SCALES_20260919.md)
reduce native-tool flags but can also conceal brief real interference.
No historical flag is removed, no cause is established and no timing is
retroactively admitted. The [completed overhead panel](DUAL_OVERHEAD_RESULT_21920.md)
has 16 validated tasks, two failed periodic measurements, and seven of nine
available pairs. The complete-panel budget remains null; the result does
not establish any comparative resource advantage.

| Original work package | Current evidence | What still proves completion |
| --- | --- | --- |
| 1. Frozen publication baseline | Comparator table, scoring corrections, OrthoMCL final-group completion and failure audit, prospective method pin | Consolidated raw-output provenance for every retained row; exact commands/versions/resources and complete claim/endpoints freeze |
| 2. Independent generalization | Completed frozen YGOB evaluation with native, command/conversion, independent-reference and overlap gates; pair enumeration and paired uncertainty | Bounded novel-taxon claim only; stronger family independence unproven; any outcome-driven method changes require new independent confirmation |
| 3. HMM and phylogeny contributions | Completed OrthoBench, original-release QfO and [corrected-QfO eight-cell factorials](QFO_CORRECTED_FACTORIAL_COMPLETE_20260919.md); paired SwissTrees simple effects/interactions; paired sequence controls; unconstrained-membership controls | Better-matched search sensitivity/calibration and controlled resource evidence remain needed. Corrected SwissTrees adjusted F1 intervals support reconciliation at C1 and C-by-R interactions, conditional on development-exposed families; all profile-refinement F1 intervals include zero. Corrected native/replay partition equality is admitted; historical-input scores cannot be transferred |
| 4. Uncertainty and error explanation | Paired OrthoBench intervals and feature strata; recovered-stage and eight-method SwissTrees intervals; corrected sequence-control intervals; annotation-defined SwissTrees domain strata and figure; native GO/EC/FAS arithmetic audits; VGNC cross-block dependency audit; TreeFam pooled count audit; initial-edge and later-stage traces; observed prefilter/scoring rejection recount and final-grouping join | Appropriate uncertainty for other QfO endpoints/secondary mean; independent duplication and fragment annotations and remaining divergence/composition analyses. Observed rejection localization does not prove causal mechanisms, historical raw-score identity, counterfactual recovery or independent FAS validation |
| 5. Robustness and practical efficiency | Corrected fixed/variable-length multi-seed simulations, complete generating-tree/NNI simulation controls, six OrthoBench tree perturbations and six parameter variants; all27 DGX repeated runs with native validation, resource replay, descriptive summaries and figure; [replacement protocol](DGX_SCALING_REPLACEMENT_PROTOCOL_20260920.md), [long-run amendment](DGX_SCALING_LONG_RUN_AMENDMENT_20260920.md) and tested measurement audit | Controlled resource comparison remains unproven. The replacement plan exists but has not run; environmental-policy freeze, authorization, launch/session integration and complete native/environmental admission remain required. QfO robustness and broader evolutionary realism remain limited |
| 6. Biological usefulness | [Prospective WGD protocol](BIOLOGICAL_WGD_APPLICATION_PROTOCOL_20260917.md), admitted runs, [results](BIOLOGICAL_WGD_RESULTS_20260917.md), [separate rescore](biological_wgd_rescore_audit_20260917.json), [figure](figures_wgd_application_20260917/wgd_application.pdf), [all six examples](figures_wgd_application_20260917/PRESPECIFIED_EXAMPLES.md), and [completed case trace](BIOLOGICAL_WGD_CASE_TRACE_20260917.md): phylogenetic OrthoHMM improves supported separation over high sensitivity but trails full OrthoFinder and SonicParanoid in separation and coverage | Bounded application and stage localization complete; tree correctness and upstream effects unresolved. Development-exposed, not independent validation or copy-specific orthology truth; YGOB column order does not identify ancestral copies; configurations differ beyond reconciliation |
| 7. Publication package | Eighteen bundled figure panels including descriptive DGX resources and corrected search controls; [corrected factorial supplement](PUBLICATION_CORRECTED_FIGURE_BUNDLE_20260919.md); [method diagram](figures_publication_method_20260916/publication_method.pdf); evolving [manuscript](PUBLICATION_MANUSCRIPT_DRAFT_20260916.md); [refreshed and relocated direct-evidence bundle](PUBLICATION_FIGURE_BUNDLE_20260918.md); bounded statistical reproduction including the corrected factorial; [current-source unit and CLI checks](PUBLICATION_TEST_REFRESH_20260918.md) | Complete corrected-QfO comparator results and manuscript integration, verified bibliography, portable executable workflows/dependencies, transitive raw-data provenance and rights clearance, versioned release and external archive. Local figure bundles and statistical exports are not full scientific reproduction |

## Current Execution Status

- The prospective replacement scaling panel is prepared, not executed. Its
  [v2 plan](dgx_root_context_scaling_plan_v2_20260920.json) retains execution
  authorization as false and the environmental policy as unresolved.
  The [non-CPU assessment](DGX_NON_CPU_ASSESSMENT_20260920.md) and positive
  native pressure observations do not justify an outside-interference claim
  or outcome-dependent exclusion. The composed measurement audit checks task
  records, raw replay and native outcomes; it neither validates outputs nor
  authorizes submissions or admits scientific timings. No service change
  has been assumed approved. Corrected BLAST 21713 remains running;
  FastOMA 21740, parameter array 21932 and CPM 21956 remain pending as of
  the September 20 scheduler check. These states are observations, not
  predicted completion dates.

- [Read-crossing control 22018](LINEAGE_READ_CROSSING_RESULT_22018.md)
  completed with all three trials and their replay passing. All signed
  partial-span values remain, including negative values. No native timing
  or historical flag was changed.

- [Lineage overhead array 21999](LINEAGE_OVERHEAD_RESULT_21999.md) and
  recorder 22000 completed. Independent replay checked 2,413 polls with zero
  observation errors; the full native audit validated all 18 tasks. All nine
  overhead pairs met numerical budgets. Interval flags remain and scientific
  timing validity is not established.

- [Dual-collector overhead array 21920](DUAL_OVERHEAD_RESULT_21920.md) and
  terminal recorder 21922 are complete. The audit retains all 18 tasks:
  16 validated, two failed, seven available pairs, and 23 narrow interval
  flags in validated periodic measurements. Median numerical budgets pass
  for high-sensitivity OrthoHMM and full OrthoFinder; satellite_v2 and the
  complete-panel budget remain unestimable. No environmental or scientific
  timing admission follows. See the [compact audit summary](dual_overhead_summary_21920_20260919.json).

- The [replacement pressure-overhead panel21889](DGX_PRESSURE_OVERHEAD_AUDIT_21889_20260919.md)
  completed with 16 validated tasks and two retained failures. Seven of nine
  pairs are available; complete-panel overhead and controlled timing remain
  unestablished. The subsequent [dual-bracket diagnostic](DUAL_NATIVE_SUBMISSION_21912_20260919.md)
  is [complete and audited](DUAL_NATIVE_RESULT_21912_20260919.md): all three
  native outputs match the prescribed prior runs, but narrow CPU flags remain
  in21 satellite_v2 intervals and one OrthoFinder interval. All original
  flags are retained. Recorder21915 completed with all three terminal records
  and no errors. This is not the27-run scientific scaling panel or timing
  admission; remaining flags require investigation.
  The [failed21869deployment](DGX_PRESSURE_OVERHEAD_FAILURE_21869.md)
  remains retained: all18tasks failed before inference under the wrong
  interpreter. The replacement uses the pinned environment interpreter and
  fresh outputs; it does not erase failures or upgrade historical timing.

- [Retained host pressure audit](DGX_PRESSURE_DIAGNOSTIC_21838.md) summarizes
  all 18 overhead tasks without changing their admission. Seventeen observed
  windows enclose native execution; task6 is partial. Small recorded memory
  stalls and nonzero CPU/I/O stalls do not identify foreign interference or
  establish isolation. No retrospective threshold or timing promotion.

- A [search routing defect](SEARCH_GPU_ROUTING_FIX_20260918.md) was reproduced
  and fixed in the current tree: CUDA availability with no eligible target
  could bypass every scoring backend.53focused tests pass. The bound native
  publication runtime is CPU-only and was reverified; broader historical
  GPU-run impact remains unaudited. Frozen executors and reported scores
  were not changed. This is not an explanation of existing search misses.
  The [OrthoBench cache witness audit](OB_GPU_BATCH_WITNESS_AUDIT_20260918.md)
  finds eligible targets in144/144directions, excluding the all-long condition
  under complete-species-pair batching, not arbitrary historical sub-batches.

- [Corrected SwissTrees sequence descriptors](CORRECTED_SWISS_SEQUENCE_STRATA_RESULT_20260918.md)
  cover all563proteins and freeze9/9entropy bins. The14recovered accessions
  and four changed old records are explicit. This is input-only preparation,
  not a corrected stratified accuracy result or validated fragment annotation.
  The [primary-strata result](CORRECTED_SWISS_STRATA_RESULT_21981.md)
  reconstructs corrected raw counts and applies the frozen27-endpoint
  bootstrap. All27endpoints were independently reproduced numerically.
  Adjusted phylogenetic-OrthoHMM versus full-OrthoFinder intervals include
  zero in both bins, and all adjusted interactions include zero. Neither
  equivalence nor a composition-specific mechanism is established. The
  [All-method and secondary descriptive displays](swiss_descriptive_strata_20260923/scores.md)
  are now complete for the seven admitted methods; OrthoMCL and empty bins
  remain explicitly missing. These displays add no inferential claims.
  A separate [identity-stratum protocol](CORRECTED_SWISS_IDENTITY_PROTOCOL_20260923.md)
  was frozen for input-only MAFFT preparation before joining to outcomes.
  [Independent verification and descriptive results](CORRECTED_SWISS_IDENTITY_RESULT_22102.md)
  are complete for all 18 families and seven admitted methods. All methods
  have lower descriptive F1 in the lower-identity bin; no subgroup significance
  or mechanism claim follows. OrthoMCL remains missing.
  Sequence identity is not calibrated evolutionary distance; independent
  fragmentation and duplication annotations remain gaps.

- Corrected DIAMOND sequence search21789 completed0:0 in02:14:43 and native
  execution admission21790 completed0:0 in00:01:30. Numeric conversion21791
  completed0:0 in01:56:01, reporting593,510,904 all-hit and321,164,891 top100
  directed rows over984,137genes/78species. Independent tuple validation21792
  completed0:0 in01:42:57 and [admitted both checkpoints](QFO_SEQUENCE_NUMERIC_ADMISSION_20260918.md);
  hit-coverage21793 completed0:0 in14:59. Its
  [label-free result](QFO_CORRECTED_SEARCH_COVERAGE_20260918.md) shows54.1825%
  of initial HMM non-self hits also occur in DIAMOND all-hits, not matched
  biological sensitivity or full-pipeline recovery. Graph-memory review21798
  completed; its [resource decision](QFO_SEQUENCE_GRAPH_RESOURCE_REVIEW_20260918.md)
  allocated32CPUs/384GiB per arm. [All-hit graph21813](QFO_SEQUENCE_ALL_HITS_EXECUTION_20260918.md)
  completed0:0 in37:34; independent admission21823 completed0:0 in10:44,
  with[retained validation](qfo_sequence_graph_admission_all_hits_21823.json).
  [Frozen pair conversion21824](QFO_SEQUENCE_ALL_HITS_PAIRS_20260918.md)
  completed0:0 in2:34 with11,300,151expected/emitted/retained clique pairs
  and zero mapping losses. Six-endpoint scoring21825 completed0:0;
  independent score admission21826 completed0:0 in3:15 after successful31:09
  native scoring. [All-hit endpoint scores are admitted](QFO_SEQUENCE_ALL_HITS_SCORE_20260918.md).
  [Top100graph21814](QFO_SEQUENCE_TOP100_EXECUTION_20260918.md)
  completed0:0 in35:24; independent graph validation21827 completed0:0 in10:49.
  Top100pair conversion21828 completed0:0 in2:36 with11,285,357mapped pairs
  and zero losses; scoring21829 completed0:0 in30:57 and admission21830
  completed0:0 in3:15. [Top100endpoint scores are admitted](QFO_SEQUENCE_TOP100_SCORE_20260918.md).
  [Self-hit semantics](QFO_SEQUENCE_SELF_HIT_REVIEW_20260918.md) retain the
  frozen cap with no self exception. Conversion completion does not establish
  matched biological sensitivity or a method advantage. Both sequence arms
  and the initial HMM control are independently admitted. The
  [completed paired SwissTrees analysis](QFO_SEQUENCE_UNCERTAINTY_RESULT_20260918.md)
  favors initial HMM precision under six-endpoint adjustment; both F1 and
  recall intervals include zero. Other-endpoint uncertainty and independent
  confirmation remain unresolved.
- Three Kingdoms historical normalized-group pair counts independently
  reproduce for all eight methods. Historical Sonic native conversion matches
  all19853 groups/288562 proteins against its retained input copies. The raw
  Danio mismatch remains in that historical run. Matched inference21795 and
  assessment21796 completed0:0; the[matched-input result](THREE_KINGDOMS_SONIC_MATCHED_RESULT_20260918.md)
  validates7272TP/48FP/80FN, F1=0.9912758996728462 and2031/2035reference
  gene coverage. This is a contemporary matched-input result, not a causal
  test of the historical mismatch or genome-wide orthology accuracy.

- Corrected Proteinortho has independently admitted six-endpoint scores and
  4,695,385 mapped native pairs in the [partial corrected table](qfo_corrected_comparison_20260919_v4/scores.md).
  Corrected Sonic inference21710 and conversion21726 completed successfully;
  15,248,739 distinct pairs survived reference mapping. Assessment21727 and
  independent admission21728 completed successfully; a fresh frozen admission
  reproduces the Sonic report byte-for-byte. Its six endpoints are admitted,
  but paired uncertainty and the remaining methods are not. Corrected HMM21706_0,
  native admission21720 and replay preparation21722 completed successfully.
  The [native result](QFO_CORRECTED_HMM_NATIVE_RESULT_20260918.md) covers984137genes
  in391908groups; its high-sensitivity point estimates are now admitted in
  the partial corrected table. The
  hit-coverage analysis21793 has completed. Checked replay21756 and its
  independent admission21757 completed successfully; all391908final groups
  match the corrected native result. [Replay evidence](QFO_CORRECTED_REPLAY_ADMITTED_20260918.md).
  Candidate preparation21758 and admission21759 are complete. The
  [complete eight-cell score export](qfo_corrected_factorial_complete_20260919/scores/scores.md)
  admits every corrected factorial cell, including the final p1_c1_r1.
  Expanded p0_c1_r1 reconciliation21760_1, independent native admission21762
  and conversion21768 are also complete, with5,977,100native pairs and zero
  mapping losses. Scoring21779 and admission21780 completed0:0; a fresh
  frozen admission reproduces the retained receipt byte-for-byte.
  [The candidate-expansion contrast](QFO_CORRECTED_EXPANDED_RECONCILIATION_SCORES_20260919.md)
  raises observed SwissTrees/VGNC/TreeFam-A F1 but lowers GO/EC/FAS.
  The [profile-refinement contrast](QFO_CORRECTED_PROFILE_RECONCILIATION_SCORES_20260919.md)
  has five slightly lower point estimates and higher EC; scoring21783 and
  admission21784 completed successfully with a byte-identical fresh recheck.
  Final admission21788 and complete-factorial uncertainty21894 completed.
  Fresh admission/count/bootstrap executions reproduce the retained JSON
  byte-for-byte. [Corrected SwissTrees intervals](QFO_CORRECTED_FACTORIAL_COMPLETE_20260919.md)
  support positive adjusted F1 interactions at both P settings, conditional
  on 18 development-exposed families. No full-method superiority is established.
  The [exploratory tree diagnostic](QFO_CORRECTED_EXPANDED_RECONCILIATION_20260919.md)
  finds38shared nontrivial rooted clades out of76in each78-species tree.
  Candidate expansion changes downstream tree estimation, so this is an
  end-to-end contrast, not a fixed-tree causal isolation or evidence that
  either inferred tree is correct. OrthoFinder21706_1 and the remaining
  corrected reconciliation/scoring chains remain unfinished.
  [Progress ledger](PUBLICATION_PROGRESS.md).
- The [recovered QfO factorial protocol](QFO_FACTORIAL_PROTOCOL_20260917.md)
  freezes eight P/C/R cells and 42 SwissTrees comparison endpoints without
  retuning. Original-release preparation, all four reconciliations, native
  admissions, conversions and eight assessment admissions are complete.
  Final assessment/admission jobs21723/21724 completed successfully; repeated
  frozen admission reproduced the entire final report. The two baseline
  cells retain validated score reuse, not independent new observations.
  R-on evaluates native inferred pairs, not RootHOG clique pairs.
  Job21725 completed the exact shared-reference count audit and paired
  bootstrap; a fresh count audit and independent arithmetic reproduced
  the results. [Results](QFO_FACTORIAL_SWISS_RESULTS_20260918.md) show no
  adjusted F1 or C-by-R interaction interval excluding zero, but consistent
  reconciliation precision gains and recall losses. Corrected-release
  inference and its separate factorial remain unfinished. Earlier failed
  submission21669 remains retained; none of this analysis uses the DGX
  timing node or establishes controlled end-to-end resource performance.
- [GO/EC arithmetic](QFO_GO_EC_ARITHMETIC_AUDIT_20260917.md) verifies all24
  retained count/mean/interval rows within serialization bounds; Darwin's
  stderr is a Student-t95% confidence half-width. [FAS arithmetic](QFO_FAS_SAMPLE_AUDIT_20260917.md)
  verifies all12 retained sample means/SEMs and exposes unequal sampling
  coverage. Neither audit supplies family-aware method-comparison intervals
  or validates underlying annotation-score correctness.

- [Annotation-defined SwissTrees strata](SWISS_DOMAIN_STRATA_RESULTS_20260917.md)
  retain all eight methods and the frozen27-endpoint analysis. Both OrthoHMM
  modes trail full OrthoFinder in F1 in both primary bins. The phylogenetic-minus-
  sensitive F1 interval is positive in the higher-type bin, but all nine
  between-bin interaction intervals include zero. Domain causality and different
  effects between strata are not established; the repeat subset is descriptive.
- [Relocated SwissTrees reproduction](SWISS_RELOCATED_REPRODUCTION_20260917.md)
  exactly reproduces statistical JSON and Markdown from a committed-source
  export in a fresh hash-pinned environment; figure generation succeeds.
  This is a bounded statistical workflow, not complete native inference,
  raw-data scoring reproduction, license clearance or the archival release.
  The [extended relocated workflow](SWISS_DOMAIN_RELOCATED_REPRODUCTION_20260917.md)
  also reproduces the domain-stratified scientific JSON and Markdown exactly
  and generates its figure from a22-file committed export; it does not
  regenerate domain annotations or establish cross-platform equivalence.
  [OrthoBench factorial relocation](ORTHOBENCH_FACTORIAL_REPRODUCTION_20260918.md)
  now also reproduces all scientific fields exactly from a committed export:
  eight cells,70families,20,000paired draws. It reuses an existing isolated
  analysis environment and does not rerun inference, conversion or official scoring.
- [Eight-method SwissTrees paired intervals](QFO_SWISS_COMPARATOR_INTERVALS_20260917.md)
  and their [three-panel figure](QFO_SWISS_COMPARATOR_FIGURE_20260917.md)
  implement the committed24-endpoint protocol. All seven comparator-minus-full-
  OrthoFinder adjusted F1 intervals are negative. Phylogenetic OrthoHMM versus
  high sensitivity has a positive point F1 difference but its adjusted interval
  includes zero; precision is positive after adjustment. These are approximate
  conditional intervals over18development-exposed families, not independent
  confirmation or evidence about the other QfO metrics.
- [VGNC reference and prediction audit](QFO_REFERENCE_MAPPING_AUDIT_20260917.md)
  reproduces exact TP/FP/FN pairs and native endpoints for all four recovered
  stages directly from prediction databases. Unscored pairs remain excluded
  under native rules; family-level uncertainty and other competitors' complete
  prediction rescoring are not established by this audit.
  The [dependency-structure audit](VGNC_DEPENDENCY_STRUCTURE_20260918.md)
  merges shared reference proteins into16,844blocks. All TP/FN are within
  blocks, but120804/15788/121468/15804FP cross blocks in the four stages,
  with only two within-block FP per stage. This prevents treating scored
  false positives as independent single-family observations without an
  explicit dependence model. Prediction-link components are diagnostics,
  not sampling units. No VGNC intervals are admitted.
- [TreeFam-A pooled count audit](QFO_TREEFAM_COUNT_AUDIT_20260917.md) reproduces
  all four recovered-stage native endpoints. Its79,320relations carry one
  pooled case label; family-level uncertainty still requires validated original
  family mapping. Ten mapped proteins have no reference relations, not missing
  predictions. Do not transfer SwissTrees resampling assumptions.
- Recovered QfO scoring21548_0..3 and admission21584 are all COMPLETED0:0.
  [Stage results](qfo_recovered_stage_summary_20260917.json) and the SwissTrees
  interval audit are complete; uncertainty for the other challenges remains open.
- DGX timing21656 is terminal for all27 tasks. Native validation21794 completed
  successfully for all27 runs; resource replay reproduces63445 observations.
  [Post-run disposition](DGX_POSTRUN_ADMISSION_20260918.md) retains all runs as
  descriptive evidence with all27 host classifications inconclusive. Process
  read errors affect4706 resource samples across18 runs. The
  [resource figure](figures_dgx_descriptive_20260918/dgx_descriptive_resources.pdf)
  shows individual values and three-repeat medians/ranges, not confidence
  intervals or controlled speedup ratios. No automatic reruns or retrospective
  monitor-rule changes are authorized by this evidence. The original32CPU
  plan remains distinct from the executed20CPU DGX plan.
  Subsequent [counter-native smokes](DGX_COUNTER_NATIVE_SMOKE_RESULT_20260918.md)
  passed native group/pair checks for all three pipelines on the645-protein
  fixture. Their counter windows cannot support host-minus-native subtraction.
  A separate [bracketed positive-control experiment](DGX_BRACKETED_CONTROLS_RESULT_20260918.md)
  passed the quiet expectation and detected both known sibling CPU loads.
  Subsequent [complete-command interval integration](DGX_INTERVAL_NATIVE_RESULT_20260918.md)
  captured all three native commands but retained adverse interval flags.
  The [residual diagnosis](DGX_INTERVAL_RESIDUAL_DIAGNOSIS_20260918.md)
  did not explain away either flag. [Hierarchy controls21816](DGX_CPU_HIERARCHY_RESULT_20260918.md)
  distinguish native sleeping-step CPU from known completed batch work with
  parent/child counter checks. These controls do not establish native-pipeline
  accounting precision, general overhead or a scientific timing inclusion
  protocol. No replacement scaling panel is authorized by those controls.
- Biological application scores, independent arithmetic, figure and all six
  prospective cases are complete. The [stage trace](BIOLOGICAL_WGD_CASE_TRACE_20260917.md)
  localizes five focal homolog separations to root-lineage grouping. Tree
  correctness and upstream search effects remain unresolved.
- [Docs dependency remediation](DOCS_DEPENDENCY_SECURITY_20260917.md) checks
  all21 retained alerts against the updated lock with zero affected versions.
  Later analysis-environment alerts were separately addressed by the
  [Pillow patch and reproduction audit](SWISS_ANALYSIS_SECURITY_20260917.md).
  The [retained read-only API snapshot](dependency_alerts_claim_audit_20260917.json),
  retrieved2026-09-18T01:38:40Z (17September local), reports zero open repository
  alerts. Earlier snapshots, including13open alerts immediately after the patch,
  remain historical. This is not a host or inference-runtime security audit
  or a current zero-alert claim: the September 20 push reports 11 dependency
  vulnerabilities (3 high, 7 moderate, 1 low). The authenticated September 20
  [snapshot and range recheck](PUBLICATION_INSTALLER_SECURITY_20260919.md)
  reconcile all 11 to the retained historical CPU-wheel lock. Its vulnerable
  installer versions remain unchanged. The separately tested patched lock's
  recorded installation versions fall outside all 11 current advisory ranges;
  this is neither repository-alert closure nor a complete current-environment,
  build-chain or exploitability audit. Do not use the historical lock for new
  installations.
  The [strict docs rebuild](docs_strict_build_validation_20260917.json)
  passes with no diagnostics; earlier14-diagnostic output remains historical.
- The [current-source regression check](PUBLICATION_REGRESSION_VERIFICATION_20260920.md)
  at `b5a1921` passed 8,314 unit tests with nine skips in 213.47 seconds.
  The overlapping opt-in native-probe/isolated-CLI suite passed 123 tests
  without skips, including all nine skipped unit identities. This does not
  constitute the full integration suite or a new frozen baseline validation.
  Regression tests do not admit pending experiments or establish
  biological validity. Earlier full runs, including `af6a8e7` (5,066 passed,
  9 skipped), remain historical evidence.
- [Retained figure integrity](PUBLICATION_FIGURE_INTEGRITY_20260918.md)
  now checks16manifests and55outputs with no byte/hash mismatch. The detached
  helper has been exported from its frozen Git revision into the
  [relocatable direct-evidence bundle](PUBLICATION_FIGURE_BUNDLE_20260918.md):
  110files,91distinct direct dependencies, verified after archive extraction.
  This does not certify scientific correctness, transitive raw data, licensing,
  full executable reproduction or external archival publication.

## Historical Execution Gates

The entries below retain prior failures and decisions. Statements describing
then-running jobs or pending analyses are historical, not current job status;
use the current status above and latest progress ledger for live work.

- YGOB `20917` was cancelled before execution after finding the missing native
  profile runtime. Corrected job `21192` completed 0:0; the frozen evaluation
  and independent pair-count crosscheck now pass. Satellite_v2 F1 is 92.233654%
  versus full OrthoFinder 92.318524%; adjusted difference interval includes zero.
  No superiority or equivalence claim. Shared-node timing remains uncontrolled.
- OrthoBench factorial preparation `21161` completed successfully. The
  [prepared manifest](orthobench_factorial_prepared_20260916.json) contains four
  candidate sets and eight planned cells; their subsequent reconciliation and
  accuracy evaluation are complete as documented immediately below.
- Reconciliation array `21248` is terminal. All four native processes succeeded
  but failed the known cwd-dependent postflight check. Independent integrity,
  native provenance and root-HOG conversion audits passed with original failed
  scheduler records retained. The complete eight-cell OrthoBench analysis and
  official-score crosschecks are in the
  [factorial results](ORTHOBENCH_FACTORIAL_RESULTS_20260916.md). Original-input
  QfO ablations are also complete; corrected-input R-off scores are admitted
  while R-on inference and scoring remain pending. OrthoBench [sequence controls](OB_SEQUENCE_SEARCH_RESULTS_20260916.md)
  and [unconstrained diagnostic](ORTHOBENCH_UNCONSTRAINED_RESULTS_20260916.md)
  are complete, without an established F1 advantage for either proposed change.
- Replay `20919` was cancelled before execution and replaced by `21088`,
  using the identical pinned command without the unnecessary YGOB dependency
  or exclusive allocation. **21088 failed equivalence**: profile construction
  silently produced zero profiles because `pair_align.so` was absent. See the
  [runtime audit](PROFILE_RUNTIME_FAILURE_20260916.md). Factorial preparation
  was blocked until the corrected-runtime replay documented below passed.
- Corrected native build is recorded in
  [runtime manifest](publication_native_runtime_20260916.json). All three CPU
  libraries built and the exact-checkout profile smoke passed. Corrected replay
  `21138` completed successfully: [all four partitions match byte-for-byte](ob_native_replay_verification_20260916.json).
  This establishes cached-stage equivalence. The OrthoBench factorial and frozen
  YGOB evaluation subsequently completed; QfO equivalence remains unresolved.
  Original defective outputs remain preserved.
- Original simulation OrthoHMM runs used the same incomplete checkout. Those results
  are defective-runtime diagnostics, not publication estimates of the intended
  method. Retain valid comparator outputs and rerun OrthoHMM with a prospectively
  recorded corrected runtime. Source hashes alone did not validate execution.
- [Corrected fixed-length results](SIMULATION_FIXED_NATIVE_RESULTS_20260916.md)
  are now assembled: 70 high-sensitivity and 64 satellite_v2 successes, six
  persistent species-tree failures, and no admitted OrthoFinder comparisons.
  Profile construction ran but added zero graph edges. This is not a superiority
  result.
- [Corrected variable-length results](SIMULATION_VARIABLE_NATIVE_RESULTS_20260916.md)
  are also assembled: 70 high-sensitivity, 67 satellite_v2 and 65 full
  OrthoFinder successes. Full OrthoFinder leads all paired condition means;
  failures and complete-case exclusions remain explicit. This does not support
  general OrthoHMM superiority. Profile expansion again added no graph edges.
- The label-blind verifier checks inference files and completion; it does
  not yet certify every scientific scoring gate. Its output retains that
  distinction explicitly.
- YGOB accuracy has now been inspected after admission. Do not change its
  frozen settings; subsequent outcome-informed changes need new confirmation.
- Dependency security alerts remain untriaged. Resolve and document them
  for release without silently changing the historical benchmark environment.
- The six prespecified OrthoBench supplied-tree perturbations completed native
  validation and official-score crosschecks. F1 ranged73.744102-74.270946%
  versus control74.106074%; all18 adjusted endpoint intervals include zero.
  [Results](OB_SPECIES_TREE_ROBUSTNESS_RESULTS_20260916.md) support small observed
  changes in this fixed exploratory panel, not equivalence or arbitrary-tree
  robustness. The OrthoBench parameter panel subsequently completed; simulation-
  truth tree-error controls and QfO robustness remain outstanding.
- QfO capture21305 confirms identical initial graph arrays/gene order but a
  different first clustering partition relative to diagnostic21295. The
  [drift diagnosis](QFO_REPLAY_DRIFT_DIAGNOSIS_20260916.md) localizes observed
  divergence before profiles, without establishing its specific cause or
  historical replay equivalence. Preserve all runs; no best-repeat selection.
- Archival deposition, submission, and any external permissions remain
  unexecuted. Document these explicitly; do not invent identifiers or approvals.
- The [84-endpoint feature-stratified analysis](OB_STRATIFIED_ERROR_RESULTS_20260916.md)
  completed with full-reference sufficient-statistic and official-score
  agreement. None of22 adjusted F1 intervals excludes zero. Two precision
  advantages and seven recall deficits survive adjustment across overlapping
  strata; neither subgroup superiority nor a causal mechanism is established.
  The one-family composition bin has no intervals; empty bins remain visible.
- QfO21307's three matching single-CPU repeats were followed by an affinity
  experiment21311 that demonstrated same-affinity partition disagreement.
  Boundary diagnostics21315/21321 found native endpoint changes before the
  optimizer. Construction-only21323 completed and was independently admitted:
  five workers matched; one explicit-int64 worker had six endpoint mismatches
  with intact original/converted arrays. Full witness/hash consistency and257
  file records were checked. This does not establish a library or hardware cause.
  [Direct-stage21326](QFO_DIRECT_GRAPH_RESULTS_20260916.md) completed and was
  independently admitted: four workers had pre-weight mismatches, including
  minimal-import workers; two matched. Differences were unchanged after weights.
  No library/hardware cause is established. [Constructor-format21327](QFO_CONSTRUCTOR_FORMAT_RESULTS_20260916.md)
  completed and passed independent admission: all three Python-pair workers match,
  one NumPy worker mismatches. This is not a proven production fix. Integrity-gated
  optimizer repeats21328 subsequently completed: all three preserve the complete
  graph before/after optimization and yield byte-identical349,898-group partitions
  covering976,504 genes. Independent admission checked308 provenance records and
  every pairwise partition comparison. This is initial-graph evidence, not a full
  cached replay or general repeatability proof. No accuracy selection has occurred.
  Full checked replay21329 later failed the profile_base worker-environment gate:
  OMP32 was inherited from profile expansion where1 was required. Its failed
  evidence is [preserved](QFO_CHECKED_REPLAY_ENVIRONMENT_FAILURE_20260917.md).
  Corrected child-only environment isolation is running as21333; no completed
  full-replay admission or QfO ablation conclusion follows yet.
- All six [OrthoBench parameter variants](OB_PARAMETER_NEIGHBORHOOD_RESULTS_20260916.md)
  passed native admission and official-score checks with251,378 genes retained.
  F1 ranges71.463468-74.973144% versus control74.106074%; all six adjusted F1
  intervals include zero. Only the CPM0.12 recall deficit excludes zero across
  the18 planned endpoints. No new default, equivalence or superiority claim.

- The [first matched scaling protocol](MATCHED_SCALING_PROTOCOL_20260916.md)
  and [input manifest](publication_scaling_inputs_20260916.json) are frozen:
  nested4/8/12 complete proteomes contain73,266/165,168/251,378 proteins, with27
  planned runs across three principal methods. No timing run has started;
  command/native-runtime, resource-accounting and host-workload gates remain
  required. This preparation does not establish practical efficiency.

- The [simulation generating-tree protocol](SIMULATION_TREE_CONTROL_PROTOCOL_20260917.md)
  and [prepared inputs](simulation_tree_controls_prepared_20260917.json) cover
  all70 variable-length datasets:210 trees and420 planned supplied-tree runs.
  Taxa/hashes/projected generating clades/distances passed independent rereading.
  Native parser testing required [plain-Newick derivatives](simulation_portable_trees_prepared_20260917.json)
  for OrthoFinder. The [unchanged-tree pilot](SIMULATION_TREE_MODE_PILOT_20260917.md)
  passed independent native/pair/topology/artifact admission for baseline_seed1
  only. Array21334 now covers all70 mode-control dataset slots:130new runs,
  2pilot runs reused and8unavailable controls from original failures. Missing
  controls remain in the future420-run oracle inventory. These are oracle and
  topology-stress inputs, not completed robustness results; full-panel mode
  admission, main supplied-tree inference and scoring remain outstanding.

- The [full checked QfO replay](QFO_CHECKED_FULL_REPLAY_RESULTS_20260917.md) now
  passed independent retrospective admission after an exact wrapper-label-failure
  recovery. Allfour clustering stages and complete976504-gene coverage passed;
  the initial partition matches checked repeats. The final390980groups differ
  from historical390817groups (4652historical-only,4815replay-only). Historical
  scores must not be transferred to this output. No full-pipeline determinism,
  accuracy benefit, or isolated runtime claim follows. Original21333FAILED1:0
  remains preserved; audit21480completed without rerunning inference.

- Superseding the earlier simulation preparation status: all132available
  unchanged-tree controls and199downstream partition comparisons passed, with
  eight unavailable originals retained. Main420method runs have finished as
  cells21405_0 and21406_1-209; independent validator21435 completed0:0. The
  [complete results](SIMULATION_TREE_ROBUSTNESS_RESULTS_20260917.md) retain405
  admitted supplied outputs and15 native failures. Cross-arm audit and all126
  exploratory endpoints are complete. Two upstream-different comparisons remain
  explicit, and no failed accuracy is imputed. All126 endpoints are now shown
  in a source-validated figure and integrated into manuscript Methods/Results;
  earlier preparation-only descriptions are historical.

- All four recovered QfO stage pair files were prepared under a frozen
  conversion/mapping protocol, with removed unmapped pairs recorded separately.
  The native assessment environment is pinned, including reference data,
  workflow, Java runtime and local container images. Scoring array21548 is
  active; independent auditor21584 waits for all four terminal stages. Neither
  successful pair conversion nor submitted native assessments establishes
  accuracy. The [stage assessment protocol](QFO_RECOVERED_STAGE_ASSESSMENT_PROTOCOL_20260917.md)
  fixes four contrasts and six endpoints; this is not the full candidate-by-
  phylogeny factorial. Appropriate paired uncertainty remains outstanding.

The manuscript must retain negative and neutral findings. Missing evidence
cannot be replaced with a claim that the method is generally applicable,
that a competitor failed, or that a custom mean proves superiority.
