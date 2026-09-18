# Publication Claim-To-Evidence Checklist

Status updated 17 September 2026. This is a completion audit, not a replacement
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
| Broad candidates improve reconciliation | [Completed OrthoBench factorial](ORTHOBENCH_FACTORIAL_INTERPRETATION_20260916.md) | Recall rises and precision falls; candidate-expansion F1 intervals include zero; QfO and additional controls pending |
| The frozen method was evaluated on novel taxa | [YGOB evaluation](YGOB_FROZEN_INTERPRETATION_20260916.md) | Complete bounded transfer evaluation; satellite F1 difference interval includes zero, precision higher and recall lower; not family-disjoint or superiority evidence |
| Validation is family-disjoint | [Homology screen](ygob_homology_screen_20260916.json) | Not established; substantial detected overlap |
| OrthoMCL final-group QfO scoring is complete | [Verified result snapshot](publication_comparison_orthomcl_complete_20260916.json) | Supported by terminal success, six native results, metadata and pair-file hashes |
| OrthoMCL BLAST failures have negligible impact | [Failure-impact audit](ORTHOMCL_FAILURE_IMPACT_20260916.md) | Not established; direct exposure is measured, indirect and counterfactual effects are not |
| Three Kingdoms demonstrates proteome-wide accuracy | [Supplementary score record](three_kingdoms_parity_20260907.json) | Unsupported; restricted BUSCO-reference universe |
| OrthoHMM is faster or more memory efficient under matched conditions | [Progress ledger](PUBLICATION_PROGRESS.md) | Not established; historical resource-accounting and workload differences remain |
| Outer PATH records prove OrthoFinder's historical companion-tool versions | [Child-PATH audit](DGX_SCALING_MIGRATION_20260917.md) | Unsupported: installed OrthoFinder rewrites its subprocess PATH; current reconstruction resolves bundled DIAMOND2.0.13/FastTree2.1.11/MCL14-137 instead of outer versions. Historical exec-path evidence still requires audit |
| ARM and x86 scoring are generally equivalent | [Portability diagnostic](native_scoring_portability_20260917.json), [development fix](NATIVE_BANDING_FIX_20260917.md), [one pipeline fixture](dgx_orthohmm_pipeline_smoke_20260917.json) | Not established: frozen narrow-band discrepancies remain; tested correction is in development source only, not the baseline; default64 matches tested synthetic fixtures and both OH modes match one simulation fixture only |
| A nearby parameter choice improves frozen-method OrthoBench F1 | [Six-variant panel](OB_PARAMETER_NEIGHBORHOOD_RESULTS_20260916.md) | Not established: all six adjusted F1 intervals include zero; no default promotion |
| The QfO native graph is reproducibly constructed | [Checked initial-graph repeats](QFO_CHECKED_REPEAT_RESULTS_20260917.md) | Three checked Python-pair runs preserve the full graph and yield identical partitions; not general determinism or complete historical replay equivalence |
| Profile-branch processing improves recovered-stage SwissTrees accuracy | [Paired SwissTrees intervals](QFO_SWISS_INTERVALS_20260917.md) | Not established: both profile contrasts have negative observed F1 differences and adjusted intervals including zero; effects occur in CASP and GH14 only. No superiority or equivalence claim |
| Recovered-stage QfO uncertainty is fully characterized | [Count audit](qfo_swiss_counts_20260917.json), [paired intervals](qfo_swiss_intervals_20260917.json) | Only SwissTrees completed: 18-family paired resampling, all12 adjusted intervals include zero. Other challenges and secondary mean require separate methods |
| Supplying the generating tree improves simulation F1 | [Complete tree panel](SIMULATION_TREE_ROBUSTNESS_RESULTS_20260917.md) | Not established: all adjusted generating-versus-inferred intervals include zero; supplied-tree completion rescues three OrthoHMM baselines but does not supply their missing inferred accuracy |
| OrthoHMM is insensitive to species-tree error | [Complete tree panel](SIMULATION_TREE_ROBUSTNESS_RESULTS_20260917.md) | Unsupported: NNI2 F1 and recall deficits have adjusted intervals below zero in four conditions; bounded exploratory result, not arbitrary-tree robustness |
| The package is publication-ready | All sections below | Not achieved |

## Completion Requirements

| Original work package | Current evidence | What still proves completion |
| --- | --- | --- |
| 1. Frozen publication baseline | Comparator table, scoring corrections, OrthoMCL final-group completion and failure audit, prospective method pin | Consolidated raw-output provenance for every retained row; exact commands/versions/resources and complete claim/endpoints freeze |
| 2. Independent generalization | Completed frozen YGOB evaluation with native, command/conversion, independent-reference and overlap gates; pair enumeration and paired uncertainty | Bounded novel-taxon claim only; stronger family independence unproven; any outcome-driven method changes require new independent confirmation |
| 3. HMM and phylogeny contributions | Corrected-runtime OrthoBench replay equivalence, completed eight-cell factorial, sequence-search replacement and unconstrained-membership controls, native coverage and incremental cost records | QfO reproducible baseline and corresponding ablations; better-matched search sensitivity/calibration; controlled resource evidence |
| 4. Uncertainty and error explanation | Paired OrthoBench intervals, individual QfO endpoints, completed feature strata, full-family stage trace and reference-incident reconciliation reconstruction | QfO strata/appropriate uncertainty; independent duplication/domain/fragment annotations; additional initial-search and rejected-edge tracing; causal explanations remain unproven |
| 5. Robustness and practical efficiency | Corrected fixed/variable-length multi-seed simulations, complete generating-tree/NNI simulation controls, six OrthoBench tree perturbations and six parameter variants, explicit failures and resource caveats | QfO robustness, matched scaling and repeated timings; broader evolutionary realism remains limited |
| 6. Biological usefulness | [Prospective WGD protocol](BIOLOGICAL_WGD_APPLICATION_PROTOCOL_20260917.md), admitted runs, [results](BIOLOGICAL_WGD_RESULTS_20260917.md), [separate native-membership rescore](biological_wgd_rescore_audit_20260917.json), [figure](figures_wgd_application_20260917/wgd_application.pdf) and [all six examples](figures_wgd_application_20260917/PRESPECIFIED_EXAMPLES.md): phylogenetic OrthoHMM improves supported separation over high sensitivity but trails full OrthoFinder and SonicParanoid in supported separation and homolog coverage | Stage-level mechanism tracing remains; development-exposed application, not independent validation or copy-specific orthology truth; YGOB column order does not identify ancestral copies; configurations differ beyond reconciliation |
| 7. Publication package | Comparison, uncertainty, ablation, simulation, strata, tree/parameter and biological application figures; [method diagram](figures_publication_method_20260916/publication_method.pdf); evolving [manuscript](PUBLICATION_MANUSCRIPT_DRAFT_20260916.md) | Scaling figures, completed Methods/Results and claim audit, verified bibliography, portable workflows/dependencies, versioned release and archival bundle |

## Current Execution Status

- The [recovered QfO factorial protocol](QFO_FACTORIAL_PROTOCOL_20260917.md)
  freezes eight P/C/R cells and 42 SwissTrees comparison endpoints without
  retuning. Candidate preparation job21670 completed successfully on bizon
  from pinned executor bd5229d; [all four arms](QFO_FACTORIAL_PREPARATION_20260917.md)
  preserve the full input universe. Reconciliation array21671 now has task0
  RUNNING and tasks1-3 pending, using executor de3202f and the verified frozen
  core launcher. No reconciliation completion or new accuracy is claimed.
  Initial submission21669 failed a mistyped commit-argument check before
  preparation and is retained. The DGX timing node is not used.
  Admission array21673 waits for reconciliation termination and requires
  native pair integrity as well as complete group coverage. A pre-scoring
  protocol clarification corrects the inherited RootHOG conversion wording:
  R-on QfO evaluates native pairwise predictions, not RootHOG clique pairs.
  No new accuracy score or inference configuration changed in that correction.
  Conversion array21674(R-off) completed all four cells;21675(R-on) waits on
  admission and independently revalidates native inferred pairs. Baseline
  assessment jobs21679/21680 completed exact independently revalidated reuse
  of recovered stages1/3, retaining original FAS samples and participant IDs.
  Expanded R-off scoring21681 is running and21682 queued. These reused scores
  are not independent new observations; fresh outcomes and the complete
  eight-cell comparison remain pending.
  Independent fresh-score admission21683/21684 waits on the corresponding
  scoring jobs and requires all native tasks, metrics and provenance to pass;
  queued validators are not evidence of completed scoring or valid intervals.
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
- [TreeFam-A pooled count audit](QFO_TREEFAM_COUNT_AUDIT_20260917.md) reproduces
  all four recovered-stage native endpoints. Its79,320relations carry one
  pooled case label; family-level uncertainty still requires validated original
  family mapping. Ten mapped proteins have no reference relations, not missing
  predictions. Do not transfer SwissTrees resampling assumptions.
- Recovered QfO scoring21548_0..3 and admission21584 are all COMPLETED0:0.
  [Stage results](qfo_recovered_stage_summary_20260917.json) and the SwissTrees
  interval audit are complete; uncertainty for the other challenges remains open.
- DGX matched timing21656 has six completed tasks and task6 running at the
  2026-09-17 accounting check. No scientific timings have yet passed output and
  resource admission. The original32CPU plan remains distinct from the DGX plan.
  [Run00 host review](DGX_RUN00_HOST_REVIEW_20260917.md) reproduces inconclusive
  observations caused by unmatched kworker-named identities. Low observed
  persistent CPU use does not certify absent contention; no monitor rule changed.
  [First-six metadata checks](DGX_FIRST_SIX_METADATA_20260917.md) match frozen
  commands, inputs and recorded runtime/resource settings, but all six retained
  host summaries are inconclusive and no native/resource admission is implied.
- Biological application scores, independent arithmetic, figure and all six
  prospective cases are complete. The [stage trace](BIOLOGICAL_WGD_CASE_TRACE_20260917.md)
  localizes five focal homolog separations to root-lineage grouping. Tree
  correctness and upstream search effects remain unresolved.
- [Docs dependency remediation](DOCS_DEPENDENCY_SECURITY_20260917.md) checks
  all21 retained alerts against the updated lock with zero affected versions.
  The [subsequent API recheck](dependency_alerts_recheck_20260917.json) reports
  zero open repository alerts. This is not a host or inference-runtime security
  audit. The [strict docs rebuild](docs_strict_build_validation_20260917.json)
  passes with no diagnostics; earlier14-diagnostic output remains historical.

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
