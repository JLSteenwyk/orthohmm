# Publication Claim-To-Evidence Checklist

Status updated 18 September 2026. This is a completion audit, not a replacement
for the original publication goal. A linked plan or passing unit test is not
evidence that an experiment completed or a biological hypothesis is true.

## Claim Boundaries

| Proposed statement | Evidence | Assessment |
| --- | --- | --- |
| Satellite_v2 has a higher observed OrthoBench aggregate F1 than full OrthoFinder | [Paired analysis](ORTHOBENCH_UNCERTAINTY_20260916.md) | Descriptively supported; interval includes zero; development-exposed |
| Satellite_v2 trades higher precision for lower recall on OrthoBench | [Paired analysis](ORTHOBENCH_UNCERTAINTY_20260916.md) | Supported within this benchmark; not selection-adjusted generalization |
| OrthoHMM outperforms full OrthoFinder overall | [Eight-method comparison](PUBLICATION_COMPARISON_ORTHOMCL_COMPLETE_20260916.md) | Not supported; endpoints and benchmark rankings differ |
| HMM expansion contributes in historical OrthoBench processing | [Historical component audit](HISTORICAL_PROFILE_ABLATION_20260916.md) | Descriptive +0.595610 F1 points; current factorial intervals include zero |
| Initial HMM search improves F1 over sequence-search replacement | [Completed sequence controls](OB_SEQUENCE_SEARCH_RESULTS_20260916.md) | Not established: observed HMM F1 higher, both adjusted difference intervals include zero; hit sensitivity/calibration unmatched |
| Broad candidates improve reconciliation | [OrthoBench factorial](ORTHOBENCH_FACTORIAL_INTERPRETATION_20260916.md), [original-release QfO factorial](QFO_FACTORIAL_SWISS_RESULTS_20260918.md) | No adjusted candidate-expansion F1 benefit established; original-QfO C-by-R interaction intervals include zero. Corrected-release QfO and additional controls remain pending |
| The frozen method was evaluated on novel taxa | [YGOB evaluation](YGOB_FROZEN_INTERPRETATION_20260916.md) | Complete bounded transfer evaluation; satellite F1 difference interval includes zero, precision higher and recall lower; not family-disjoint or superiority evidence |
| Validation is family-disjoint | [Homology screen](ygob_homology_screen_20260916.json) | Not established; substantial detected overlap |
| Original-release OrthoMCL final-group QfO scoring is complete | [Verified result snapshot](publication_comparison_orthomcl_complete_20260916.json) | Supported for the original inputs only; corrected-release BLAST is queued and no corrected OrthoMCL score is available |
| OrthoMCL BLAST failures have negligible impact | [Failure-impact audit](ORTHOMCL_FAILURE_IMPACT_20260916.md) | Not established; direct exposure is measured, indirect and counterfactual effects are not |
| Three Kingdoms demonstrates proteome-wide accuracy | [Supplementary score record](three_kingdoms_parity_20260907.json) | Unsupported; restricted BUSCO-reference universe |
| Every historical Three Kingdoms method used identical input bytes | [Method-input audit](THREE_KINGDOMS_METHOD_INPUTS_20260918.md), [matched rerun](THREE_KINGDOMS_MATCHED_SONIC_20260918.md) | Not established: SonicParanoid native snapshot matches raw Danio rather than the staged stop-marker-stripped version; older high-sensitivity record lacks per-file hashes. Matched inference21795 and assessment21796 are queued, not completed |
| Historical Three Kingdoms scores reproduce from normalized groups | [Independent arithmetic audit](THREE_KINGDOMS_PAIR_COUNT_AUDIT_20260918.md) | Supported for all eight retained methods; native Sonic group conversion also verified. This does not establish matched historical inputs or proteome-wide accuracy |
| OrthoHMM is faster or more memory efficient under matched conditions | [DGX disposition](DGX_POSTRUN_ADMISSION_20260918.md), [descriptive observations](dgx_descriptive_resources_20260918.json) | Not established; all27 native-valid runs retained as descriptive evidence, with host-isolation uncertainty and asymmetric process-sampling gaps. No controlled comparison admitted |
| Outer PATH records prove OrthoFinder's historical companion-tool versions | [Child-PATH audit](DGX_SCALING_MIGRATION_20260917.md) | Unsupported: installed OrthoFinder rewrites its subprocess PATH; current reconstruction resolves bundled DIAMOND2.0.13/FastTree2.1.11/MCL14-137 instead of outer versions. Historical exec-path evidence still requires audit |
| ARM and x86 scoring are generally equivalent | [Portability diagnostic](native_scoring_portability_20260917.json), [development fix](NATIVE_BANDING_FIX_20260917.md), [one pipeline fixture](dgx_orthohmm_pipeline_smoke_20260917.json) | Not established: frozen narrow-band discrepancies remain; tested correction is in development source only, not the baseline; default64 matches tested synthetic fixtures and both OH modes match one simulation fixture only |
| A nearby parameter choice improves frozen-method OrthoBench F1 | [Six-variant panel](OB_PARAMETER_NEIGHBORHOOD_RESULTS_20260916.md) | Not established: all six adjusted F1 intervals include zero; no default promotion |
| The QfO native graph is reproducibly constructed | [Checked initial-graph repeats](QFO_CHECKED_REPEAT_RESULTS_20260917.md) | Three checked Python-pair runs preserve the full graph and yield identical partitions; not general determinism or complete historical replay equivalence |
| Profile-branch processing improves recovered-stage SwissTrees accuracy | [Paired SwissTrees intervals](QFO_SWISS_INTERVALS_20260917.md) | Not established: both profile contrasts have negative observed F1 differences and adjusted intervals including zero; effects occur in CASP and GH14 only. No superiority or equivalence claim |
| Recovered-stage QfO uncertainty is fully characterized | [Count audit](qfo_swiss_counts_20260917.json), [paired intervals](qfo_swiss_intervals_20260917.json) | Only SwissTrees completed: 18-family paired resampling, all12 adjusted intervals include zero. Other challenges and secondary mean require separate methods |
| OrthoHMM phylogeny exceeds full OrthoFinder on SwissTrees F1 | [Eight-method paired intervals](QFO_SWISS_COMPARATOR_INTERVALS_20260917.md) | Unsupported: its adjusted F1 difference interval is negative; conditional evidence from 18 development-exposed families, not all QfO endpoints |
| Domain architecture explains the OrthoHMM configuration effect | [Stratified results](SWISS_DOMAIN_STRATA_RESULTS_20260917.md), [all27 endpoints plotted](SWISS_DOMAIN_STRATA_FIGURE_20260917.md) | Not established: all nine adjusted interaction intervals include zero; annotation-defined association does not establish causality |
| Satellite constraints explain the five focal WGD homolog-coverage losses | [Prespecified case trace](BIOLOGICAL_WGD_CASE_TRACE_20260917.md) | Contradicted for these five homologs: all were in anchor candidates and separated during root-lineage reconstruction before constraints; topology correctness and upstream effects remain unresolved |
| The original-release QfO factorial has paired uncertainty estimates | [Validated results](QFO_FACTORIAL_SWISS_RESULTS_20260918.md), [counts](qfo_factorial_swiss_counts_20260918.json), [42 endpoints](qfo_factorial_swiss_bootstrap_20260918.json) | Complete for SwissTrees only: all adjusted F1 intervals include zero; R increases precision and lowers recall. No adjusted C-by-R interaction established. Corrected-release reruns remain separate and unfinished |
| Native GO/EC/FAS error bars establish paired method differences | [GO/EC audit](QFO_GO_EC_ARITHMETIC_AUDIT_20260917.md), [FAS sample audit](QFO_FAS_SAMPLE_AUDIT_20260917.md) | Unsupported: GO/EC use Student-t95% half-widths, FAS uses sample SEM, and none supplies dependency-aware paired method intervals |
| Supplying the generating tree improves simulation F1 | [Complete tree panel](SIMULATION_TREE_ROBUSTNESS_RESULTS_20260917.md) | Not established: all adjusted generating-versus-inferred intervals include zero; supplied-tree completion rescues three OrthoHMM baselines but does not supply their missing inferred accuracy |
| OrthoHMM is insensitive to species-tree error | [Complete tree panel](SIMULATION_TREE_ROBUSTNESS_RESULTS_20260917.md) | Unsupported: NNI2 F1 and recall deficits have adjusted intervals below zero in four conditions; bounded exploratory result, not arbitrary-tree robustness |
| The package is publication-ready | All sections below | Not achieved |
| Original QfO inputs match the corrected2020benchmark release | [Corrected archive comparison](QFO_CORRECTED_ARCHIVE_ACQUIRED_20260918.md) | Contradicted for the Xenopus proteome; preserve original results as release-limited |
| Corrected QfO inputs cover the retained reference identities and sequence content | [Native sequence and staging audits](QFO_CORRECTED_INPUTS_STAGED_20260918.md) | Supported: all 984,137 identities, 983,959 exact sequences and 178 representation-only differences; no unexplained differences. This is input compatibility, not biological annotation validation |
| Corrected QfO accuracy or rankings are established | [Admitted partial table](qfo_corrected_comparison_20260918_v2/scores.md), [frozen rerun protocol](QFO_CORRECTED_RELEASE_PROTOCOL_20260917.md) | Proteinortho and SonicParanoid point estimates are admitted; the eight-method comparison and paired differences remain incomplete. Remaining methods and factorial are unfinished. No original predictions or intervals may be relabeled |
| Original TreeFam family-level uncertainty can be recovered from pooled pairs | [Source retrieval investigation](TREEFAM_SOURCE_RETRIEVAL_20260918.md) | Unsupported: original trees and mapping remain missing; downloaded pooled reference is not an independent-family inventory |
| Merging overlapping VGNC labels makes ordinary family resampling valid | [Dependency audit](VGNC_DEPENDENCY_STRUCTURE_20260918.md) | Not established: 16,863 labels form 16,844 reference blocks, but almost all scored false positives cross blocks. Outcome-defined prediction components are not independent reference units |
| Missing initial OrthoBench edges explain final grouping errors | [Initial-edge trace](OB_INITIAL_EDGE_TRACE_20260918.md) | Descriptive localization only: 505 hit-supported separated memberships had initial edges and 1,466 did not. Other graph paths and later stages prevent a causal conclusion |
| OrthoBench factorial statistics reproduce outside the checkout | [Relocated reproduction](ORTHOBENCH_FACTORIAL_REPRODUCTION_20260918.md) | Exact agreement for the eight-cell, 70-family, 20,000-draw statistical analysis; not native inference/scoring or cross-platform reproduction |
| New counter controls establish controlled comparative timing | [Complete native-command integration](DGX_INTERVAL_NATIVE_RESULT_20260918.md), [residual diagnosis](DGX_INTERVAL_RESIDUAL_DIAGNOSIS_20260918.md), [hierarchy controls](DGX_CPU_HIERARCHY_RESULT_20260918.md) | Not established: satellite_v2 and OrthoFinder each retain an interval flag. Read-window increments and observed leaf CPU do not explain them; separate sleeping-step controls localize known batch loads but do not resolve native-run accounting/overhead. Non-CPU isolation and scientific inclusion rules remain incomplete. No historical timing is upgraded |

## Completion Requirements

| Original work package | Current evidence | What still proves completion |
| --- | --- | --- |
| 1. Frozen publication baseline | Comparator table, scoring corrections, OrthoMCL final-group completion and failure audit, prospective method pin | Consolidated raw-output provenance for every retained row; exact commands/versions/resources and complete claim/endpoints freeze |
| 2. Independent generalization | Completed frozen YGOB evaluation with native, command/conversion, independent-reference and overlap gates; pair enumeration and paired uncertainty | Bounded novel-taxon claim only; stronger family independence unproven; any outcome-driven method changes require new independent confirmation |
| 3. HMM and phylogeny contributions | Completed OrthoBench and original-release QfO eight-cell factorials; paired SwissTrees simple effects/interactions; sequence-search replacement and unconstrained-membership controls | Corrected-release QfO factorial; better-matched search sensitivity/calibration; controlled resource evidence. Checked QfO replay differs from historical final groups; do not transfer historical scores |
| 4. Uncertainty and error explanation | Paired OrthoBench intervals and feature strata; recovered-stage and eight-method SwissTrees intervals; annotation-defined SwissTrees domain strata and figure; native GO/EC/FAS arithmetic audits; VGNC prediction rescore and cross-block dependency audit; TreeFam pooled count audit; initial-edge and later-stage traces | Appropriate uncertainty for other QfO endpoints/secondary mean; independent duplication and fragment annotations and remaining divergence/composition analyses; prefilter-versus-scoring rejection remains unresolved. Completed edge tracing and domain inventory do not prove causal mechanisms or independent FAS validation |
| 5. Robustness and practical efficiency | Corrected fixed/variable-length multi-seed simulations, complete generating-tree/NNI simulation controls, six OrthoBench tree perturbations and six parameter variants; all27 DGX repeated runs with native validation, resource replay, descriptive summaries and figure | Controlled resource comparison remains unproven despite completed repeated timing; a new prospective observation/inclusion plan is required before any controlled-speed experiment. QfO robustness and broader evolutionary realism remain limited |
| 6. Biological usefulness | [Prospective WGD protocol](BIOLOGICAL_WGD_APPLICATION_PROTOCOL_20260917.md), admitted runs, [results](BIOLOGICAL_WGD_RESULTS_20260917.md), [separate rescore](biological_wgd_rescore_audit_20260917.json), [figure](figures_wgd_application_20260917/wgd_application.pdf), [all six examples](figures_wgd_application_20260917/PRESPECIFIED_EXAMPLES.md), and [completed case trace](BIOLOGICAL_WGD_CASE_TRACE_20260917.md): phylogenetic OrthoHMM improves supported separation over high sensitivity but trails full OrthoFinder and SonicParanoid in separation and coverage | Bounded application and stage localization complete; tree correctness and upstream effects unresolved. Development-exposed, not independent validation or copy-specific orthology truth; YGOB column order does not identify ancestral copies; configurations differ beyond reconciliation |
| 7. Publication package | Sixteen retained figure panels including descriptive DGX resources; [method diagram](figures_publication_method_20260916/publication_method.pdf); evolving [manuscript](PUBLICATION_MANUSCRIPT_DRAFT_20260916.md); [relocated direct-evidence bundle](PUBLICATION_FIGURE_BUNDLE_20260918.md); bounded SwissTrees and OrthoBench statistical reproduction | Completed corrected-QfO results and manuscript integration, verified bibliography, portable executable workflows/dependencies, transitive raw-data provenance and rights clearance, versioned release and external archive. The local figure bundle and statistical exports are not full scientific reproduction |

## Current Execution Status

- [Corrected SwissTrees sequence descriptors](CORRECTED_SWISS_SEQUENCE_STRATA_RESULT_20260918.md)
  cover all563proteins and freeze9/9entropy bins. The14recovered accessions
  and four changed old records are explicit. This is input-only preparation,
  not a corrected stratified accuracy result or validated fragment annotation.
  The [primary-strata runner](CORRECTED_SWISS_STRATA_EXECUTION_20260918.md)
  now reconstructs corrected raw counts and applies the frozen27-endpoint
  bootstrap. It has passed synthetic and input-binding tests, but actual
  execution awaits admitted corrected predictions. It does not make the
  all-method or secondary-stratum displays complete.

- Corrected DIAMOND sequence search21789 completed0:0 in02:14:43 and native
  execution admission21790 completed0:0 in00:01:30. Numeric conversion21791
  completed0:0 in01:56:01, reporting593,510,904 all-hit and321,164,891 top100
  directed rows over984,137genes/78species. Independent tuple validation21792
  completed0:0 in01:42:57 and [admitted both checkpoints](QFO_SEQUENCE_NUMERIC_ADMISSION_20260918.md);
  hit-coverage21793 still waits for HMM admission. Graph-memory review21798
  completed; its [resource decision](QFO_SEQUENCE_GRAPH_RESOURCE_REVIEW_20260918.md)
  allocated32CPUs/384GiB per arm. [All-hit graph21813](QFO_SEQUENCE_ALL_HITS_EXECUTION_20260918.md)
  completed0:0 in37:34; independent admission21823 completed0:0 in10:44,
  with[retained validation](qfo_sequence_graph_admission_all_hits_21823.json).
  Frozen pair conversion21824 is running; top100
  graph21814 is running after it. Execution-reported partitions are not yet
  independently admitted or scored.
  [Self-hit semantics](QFO_SEQUENCE_SELF_HIT_REVIEW_20260918.md) retain the
  frozen cap with no self exception. Conversion completion does not establish
  matched biological sensitivity, graph feasibility or
  an accuracy result. No independently admitted sequence graph/scoring result
  is available yet.
- Three Kingdoms historical normalized-group pair counts independently
  reproduce for all eight methods. Historical Sonic native conversion matches
  all19853 groups/288562 proteins against its retained input copies. The raw
  Danio mismatch remains; matched inference21795 is running and assessment21796
  remains dependency-pending.

- Corrected Proteinortho has independently admitted six-endpoint scores and
  4,695,385 mapped native pairs in the [partial corrected table](qfo_corrected_comparison_20260918_v2/scores.md).
  Corrected Sonic inference21710 and conversion21726 completed successfully;
  15,248,739 distinct pairs survived reference mapping. Assessment21727 and
  independent admission21728 completed successfully; a fresh frozen admission
  reproduces the Sonic report byte-for-byte. Its six endpoints are admitted,
  but paired uncertainty and the remaining methods are not. Corrected HMM21706_0 is running; its native
  admission21720 and replay preparation21722 are pending. The corrected
  replay/candidate/reconciliation/conversion/scoring/admission workflows are
  implemented and tested, not completed experiments. [Progress ledger](PUBLICATION_PROGRESS.md).
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
  remain historical. This is not a host or inference-runtime security audit.
  The [strict docs rebuild](docs_strict_build_validation_20260917.json)
  passes with no diagnostics; earlier14-diagnostic output remains historical.
- The latest full unit-suite run ataf6a8e7 passed5066tests with9skipped in89.12s,
  recorded in the [progress ledger](PUBLICATION_PROGRESS.md).
  Regression tests do not admit pending experiments or establish biological
  validity. The prior5c7c19d and4e2e93c runs remain historical evidence.
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
  [factorial results](ORTHOBENCH_FACTORIAL_RESULTS_20260916.md); QfO ablations,
  remain pending. OrthoBench [sequence controls](OB_SEQUENCE_SEARCH_RESULTS_20260916.md)
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
