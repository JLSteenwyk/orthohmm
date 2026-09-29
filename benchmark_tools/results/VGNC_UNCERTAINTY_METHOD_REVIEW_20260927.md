# VGNC Uncertainty: Replacement Method Review

The [corrected all-method deletion diagnostic](CORRECTED_VGNC_INFLUENCE_RESULT_20260928.md)
finds no single reference-block deletion reverses any comparator's F1 deficit
against full OrthoFinder. This fixed-table sensitivity result supplies neither
a sampling law nor confidence intervals; the requirements below remain open.

The subsequent [corrected eight-method component audit](CORRECTED_VGNC_COMPONENTS_20260928.md)
describes all individual, paired-union and all-method-union cross-block graphs.
Largest components range from four blocks for the phylogenetic
OrthoHMM/full-OrthoFinder union to 370 for the all-method union. This
method-dependent grouping is not an independently justified sampling law;
the statistical requirements below remain unresolved.

Primary-source review on 27 September 2026. This is a methodological decision
record, not a completed estimator, statistical proof or native confidence interval.

## Existing Negative Evidence

The [rare-error screen](SPARSE_DYADIC_F1_RESULT_20260926.md) retained only
61.3-65.1% coverage in its two boundary cells, with 348-387 invalid variances
per 1,000 replicates. When the extra error count is zero, paired observed
contributions coincide even though the generating target differs. The
[earlier dependence screen](DYADIC_F1_VARIANCE_RESULT_20260926.md) also failed
under shared-clade dependence. These are separate problems: changing sparse
storage, increasing family count or clipping negative variances does not
establish valid inference for either one.

## Primary Sources And Applicability

Matsushita and Otsu, *Empirical Likelihood for Network Data*, JASA 119 (2024),
2117-2128, develop modified jackknife empirical likelihood for network moments
and address failures of the unmodified statistic under sparsity. This is a
candidate methodological direction, not a theorem for our paired VGNC F1.
The publisher abstract and indexed article text were accessible; direct
full-page retrieval failed. No claim of complete proof verification is made.
[Publisher article](https://www.tandfonline.com/doi/abs/10.1080/01621459.2023.2250091).

Hounyo, Lin and Song, *Robust Inference for Dyadic Data with Dependent Ordered
Nodes*, arXiv:2605.28349v2 (2 June 2026), study dependence between disjoint
dyads through ordered nodes. Section 2.2 treats a meaningful order as given;
Assumption 1 specifies stationary weakly dependent node variables and
independent dyad shocks. The authors explicitly leave estimated-order
uncertainty for future work. The full HTML was inspected. This is a preprint,
not a validated biological dependence model.
[Versioned full text](https://arxiv.org/html/2605.28349v2).

**Project inference:** neither source establishes validity for the retained
VGNC data. The endpoint combines within-family truth counts and between-family
false positives, then forms two nonlinear F1 ratios. We have neither a derived
mapping to the network-moment theorem nor an independently justified ordered
family dependence structure. Alphabetical family IDs or ordering families by
observed errors would not supply that justification.

## Requirements Before Another Native Interval

1. State the population target and sampling mechanism. Distinguish the
   observed fixed-reference score from inference to new families or clades.
2. Specify the joint sampling law for paired TP/FP/FN counts, including both
   within-family and cross-family terms, shared methods and unequal sizes.
3. Derive the estimator for the actual ratio difference. Establish how its
   assumptions cover sparsity, zero extra errors and any family dependence;
   a citation to a different statistic is insufficient.
4. Freeze a new validation protocol before new results. Retain the existing
   regular, unequal-size, rare-error and shared-clade cells; failures must not
   be omitted or treated as zero-width successes. Evaluate coverage, interval
   width and invalid outcomes against a known generating target.
5. Independently justify native resampling units and any dependence ordering.
   Passing simulations under an assumed model does not prove this condition.

Decision: do not replace the failed Wald candidate with a named bootstrap or
empirical-likelihood routine solely to obtain intervals. No new native CI,
test result, endpoint, default or benchmark execution is admitted here. GO/EC,
FAS and the secondary aggregate remain separate uncertainty questions.

Retrieval limits: the Wisconsin `Paper2.pdf` URL returned 404; the older LSE
`JEL12_accepted.pdf` route was unavailable through the browser. Those failures
are not evidence against the methods. The two sources above support this
bounded review; no email or external issue was sent.
