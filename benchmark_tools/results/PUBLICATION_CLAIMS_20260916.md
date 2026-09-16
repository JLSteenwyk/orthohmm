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
| HMM expansion contributes in historical OrthoBench processing | [Historical component audit](HISTORICAL_PROFILE_ABLATION_20260916.md) | Descriptive +0.595610 F1 points after refinement; no matched HMM-free control |
| Broad candidates improve reconciliation | [Completed OrthoBench factorial](ORTHOBENCH_FACTORIAL_INTERPRETATION_20260916.md) | Recall rises and precision falls; candidate-expansion F1 intervals include zero; QfO and additional controls pending |
| The frozen method was evaluated on novel taxa | [YGOB evaluation](YGOB_FROZEN_INTERPRETATION_20260916.md) | Complete bounded transfer evaluation; satellite F1 difference interval includes zero, precision higher and recall lower; not family-disjoint or superiority evidence |
| Validation is family-disjoint | [Homology screen](ygob_homology_screen_20260916.json) | Not established; substantial detected overlap |
| OrthoMCL final-group QfO scoring is complete | [Verified result snapshot](publication_comparison_orthomcl_complete_20260916.json) | Supported by terminal success, six native results, metadata and pair-file hashes |
| OrthoMCL BLAST failures have negligible impact | [Failure-impact audit](ORTHOMCL_FAILURE_IMPACT_20260916.md) | Not established; direct exposure is measured, indirect and counterfactual effects are not |
| Three Kingdoms demonstrates proteome-wide accuracy | [Supplementary score record](three_kingdoms_parity_20260907.json) | Unsupported; restricted BUSCO-reference universe |
| OrthoHMM is faster or more memory efficient under matched conditions | [Progress ledger](PUBLICATION_PROGRESS.md) | Not established; historical resource-accounting and workload differences remain |
| The package is publication-ready | All sections below | Not achieved |

## Completion Requirements

| Original work package | Current evidence | What still proves completion |
| --- | --- | --- |
| 1. Frozen publication baseline | Comparator table, scoring corrections, OrthoMCL final-group completion and failure audit, prospective method pin | Consolidated raw-output provenance for every retained row; exact commands/versions/resources and complete claim/endpoints freeze |
| 2. Independent generalization | Completed frozen YGOB evaluation with native, command/conversion, independent-reference and overlap gates; pair enumeration and paired uncertainty | Bounded novel-taxon claim only; stronger family independence unproven; any outcome-driven method changes require new independent confirmation |
| 3. HMM and phylogeny contributions | Historical four partitions; prospective eight-cell protocol; failed replay exposed missing profile runtime | Corrected-runtime replay equivalence, executable frozen factorial, all cells on OrthoBench/QfO, membership-filter diagnostic, matched sequence control, coverage and cost measurements |
| 4. Uncertainty and error explanation | OrthoBench paired intervals, family summaries, individual QfO endpoints | Prespecified label-independent strata, stage-level error tracing, remaining appropriate uncertainty and prediction coverage |
| 5. Robustness and practical efficiency | Historical timing records and known accounting caveats | Validated evolutionary simulations, multiple seeds, duplication/loss/divergence/missingness/sampling conditions, tree error and parameter neighborhood, matched scaling and repeated timings |
| 6. Biological usefulness | No completed prespecified application | Independently supported family selection, difficult positives and negatives, relevant comparators, successes and failures without outcome-based selection |
| 7. Publication package | Generated comparison/uncertainty figures, initial [manuscript draft](PUBLICATION_MANUSCRIPT_DRAFT_20260916.md), this checklist | Method diagram, ablation/error/scaling/application figures, completed Methods/Results, verified bibliography, portable workflows/dependencies, versioned release and archival bundle |

## Current Execution Gates

- YGOB `20917` was cancelled before execution after finding the missing native
  profile runtime. Corrected job `21192` completed 0:0; the frozen evaluation
  and independent pair-count crosscheck now pass. Satellite_v2 F1 is 92.233654%
  versus full OrthoFinder 92.318524%; adjusted difference interval includes zero.
  No superiority or equivalence claim. Shared-node timing remains uncontrolled.
- OrthoBench factorial preparation `21161` completed successfully. The
  [prepared manifest](orthobench_factorial_prepared_20260916.json) contains four
  candidate sets and eight planned cells; reconciliation and accuracy evaluation
  are not yet complete.
- Reconciliation array `21248` is terminal. All four native processes succeeded
  but failed the known cwd-dependent postflight check. Independent integrity,
  native provenance and root-HOG conversion audits passed with original failed
  scheduler records retained. The complete eight-cell OrthoBench analysis and
  official-score crosschecks are in the
  [factorial results](ORTHOBENCH_FACTORIAL_RESULTS_20260916.md); QfO ablations,
  matched sequence control and unconstrained membership diagnostic remain pending.
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
  This establishes cached-stage equivalence; factorial and corrected validation
  remain incomplete. Original defective outputs remain preserved.
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
- Archival deposition, submission, and any external permissions remain
  unexecuted. Document these explicitly; do not invent identifiers or approvals.

The manuscript must retain negative and neutral findings. Missing evidence
cannot be replaced with a claim that the method is generally applicable,
that a competitor failed, or that a custom mean proves superiority.
