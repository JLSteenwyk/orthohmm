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
| Broad candidates improve reconciliation | [Prospective factorial](PUBLICATION_ABLATION_PROTOCOL_20260916.md) | Untested by completed matched factorial; a plan is not a result |
| The method transfers to novel taxa | [Frozen YGOB protocol](YGOB_VALIDATION_PROTOCOL_20260916.md) | Pending successful inference and gated evaluation |
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
| 2. Independent generalization | YGOB input/reference freeze, overlap screen, source-level resource audit, queued inference, tested scorer/report helpers | Successful run, exact command/version and native conversion checks, reference reconstruction, outcomes and paired uncertainty; explicit overlap-limited claim |
| 3. HMM and phylogeny contributions | Historical four partitions; prospective eight-cell protocol; queued current-source replay check | Replay equivalence, executable frozen factorial, all cells on OrthoBench/QfO, membership-filter diagnostic, matched sequence control, coverage and cost measurements |
| 4. Uncertainty and error explanation | OrthoBench paired intervals, family summaries, individual QfO endpoints | Prespecified label-independent strata, stage-level error tracing, remaining appropriate uncertainty and prediction coverage |
| 5. Robustness and practical efficiency | Historical timing records and known accounting caveats | Validated evolutionary simulations, multiple seeds, duplication/loss/divergence/missingness/sampling conditions, tree error and parameter neighborhood, matched scaling and repeated timings |
| 6. Biological usefulness | No completed prespecified application | Independently supported family selection, difficult positives and negatives, relevant comparators, successes and failures without outcome-based selection |
| 7. Publication package | Generated comparison/uncertainty figures, initial [manuscript draft](PUBLICATION_MANUSCRIPT_DRAFT_20260916.md), this checklist | Method diagram, ablation/error/scaling/application figures, completed Methods/Results, verified bibliography, portable workflows/dependencies, versioned release and archival bundle |

## Current Execution Gates

- YGOB `20917`: pending resources at the last live check. Do not infer
  completion from a queued job or from partial output files.
- Replay `20919`: depends on successful YGOB completion. A verified replay
  is a prerequisite to using current-source cached stages in the factorial.
- The label-blind verifier checks inference files and completion; it does
  not yet certify every scientific scoring gate. Its output retains that
  distinction explicitly.
- No YGOB accuracy has been inspected. Do not change its frozen settings
  after outcomes; subsequent outcome-informed changes need new confirmation.
- Dependency security alerts remain untriaged. Resolve and document them
  for release without silently changing the historical benchmark environment.
- Archival deposition, submission, and any external permissions remain
  unexecuted. Document these explicitly; do not invent identifiers or approvals.

The manuscript must retain negative and neutral findings. Missing evidence
cannot be replaced with a claim that the method is generally applicable,
that a competitor failed, or that a custom mean proves superiority.
