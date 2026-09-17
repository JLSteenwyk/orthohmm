# Publication Claim-To-Evidence Checklist

Status on 16 September 2026. This is a completion audit, not a replacement
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
| A nearby parameter choice improves frozen-method OrthoBench F1 | [Six-variant panel](OB_PARAMETER_NEIGHBORHOOD_RESULTS_20260916.md) | Not established: all six adjusted F1 intervals include zero; no default promotion |
| The QfO native graph is reproducibly constructed | [Direct-stage diagnostic](QFO_DIRECT_GRAPH_RESULTS_20260916.md) | Not established: pre-weight mismatches occur with intact constructor arrays in both import modes; constructor-format comparison pending |
| The package is publication-ready | All sections below | Not achieved |

## Completion Requirements

| Original work package | Current evidence | What still proves completion |
| --- | --- | --- |
| 1. Frozen publication baseline | Comparator table, scoring corrections, OrthoMCL final-group completion and failure audit, prospective method pin | Consolidated raw-output provenance for every retained row; exact commands/versions/resources and complete claim/endpoints freeze |
| 2. Independent generalization | Completed frozen YGOB evaluation with native, command/conversion, independent-reference and overlap gates; pair enumeration and paired uncertainty | Bounded novel-taxon claim only; stronger family independence unproven; any outcome-driven method changes require new independent confirmation |
| 3. HMM and phylogeny contributions | Corrected-runtime OrthoBench replay equivalence, completed eight-cell factorial, sequence-search replacement and unconstrained-membership controls, native coverage and incremental cost records | QfO reproducible baseline and corresponding ablations; better-matched search sensitivity/calibration; controlled resource evidence |
| 4. Uncertainty and error explanation | Paired OrthoBench intervals, individual QfO endpoints, completed feature strata, full-family stage trace and reference-incident reconciliation reconstruction | QfO strata/appropriate uncertainty; independent duplication/domain/fragment annotations; additional initial-search and rejected-edge tracing; causal explanations remain unproven |
| 5. Robustness and practical efficiency | Corrected fixed/variable-length multi-seed simulations, six OrthoBench tree perturbations and six parameter variants, explicit failures and resource caveats | QfO robustness, simulation-truth tree-error controls, matched scaling and repeated timings; broader evolutionary realism remains limited |
| 6. Biological usefulness | No completed prespecified application | Independently supported family selection, difficult positives and negatives, relevant comparators, successes and failures without outcome-based selection |
| 7. Publication package | Comparison, uncertainty, ablation, simulation, strata, tree/parameter figures; [method diagram](figures_publication_method_20260916/publication_method.pdf); evolving [manuscript](PUBLICATION_MANUSCRIPT_DRAFT_20260916.md) | Scaling/application figures, completed Methods/Results and claim audit, verified bibliography, portable workflows/dependencies, versioned release and archival bundle |

## Current Execution Gates

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
  No library/hardware cause is established. Constructor-format21327 remains
  running; no accuracy scores or partition selection are involved.
- All six [OrthoBench parameter variants](OB_PARAMETER_NEIGHBORHOOD_RESULTS_20260916.md)
  passed native admission and official-score checks with251,378 genes retained.
  F1 ranges71.463468-74.973144% versus control74.106074%; all six adjusted F1
  intervals include zero. Only the CPM0.12 recall deficit excludes zero across
  the18 planned endpoints. No new default, equivalence or superiority claim.

The manuscript must retain negative and neutral findings. Missing evidence
cannot be replaced with a claim that the method is generally applicable,
that a competitor failed, or that a custom mean proves superiority.
