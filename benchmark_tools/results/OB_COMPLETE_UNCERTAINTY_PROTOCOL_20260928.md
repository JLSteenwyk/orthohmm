# Complete OrthoBench Comparison Protocol

Commit and push before executing this new bootstrap. Aggregate scores and
descriptive strata are already known: this is an exploratory extension, not
preregistration or independent validation. Preserve the historical three-method
intervals and all scientific settings unchanged.

Use all eight retained methods in the pinned original-three-method report
and five-comparator readback. Retain every one of the 70 RefOGs, identical
family sizes, full-reference TP/FP/FN and low-certainty conventions. Recheck
source and prediction bytes. Do not rerun inference or select method outputs.

Compare each of the seven alternatives against full OrthoFinder 3.1.5, for
weighted F1, precision and recall: 21 endpoints. Draw 100,000 shared RefOG
multinomial resamples with NumPy PCG64 seed 20260928. Recompute sums of
TP/(n-1), FP/(n-1), FN/(n-1) and then the actual statistic in every draw.
Report differences in percentage points, nominal 95% percentile intervals,
and Bonferroni percentile intervals across all 21 endpoints. Adjusted tails
contain about 119 draws; percentile coverage remains approximate.

Report all point estimates, intervals and descriptive per-family F1
wins/ties/losses, not selected winning endpoints or bootstrap p-values.
No pair-IID resampling, tuning or new subgroup analysis. Paired resampling
does not cure development exposure, dependence between families, shared
evolutionary history or cross-family prediction errors. Interpret intervals
conditionally on exchangeable reference families, not as proven population
coverage or evidence of universal superiority.

Independently reproduce the complete 21-endpoint output with an alternative
count-summing path using the same frozen draws, and check point estimates
against the retained overall scores. Numerical agreement validates arithmetic,
not the exchangeability assumption. Retain negative and neutral findings.
