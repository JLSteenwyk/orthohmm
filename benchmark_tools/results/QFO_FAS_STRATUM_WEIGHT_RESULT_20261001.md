# Retained FAS Stratum Weight Decomposition

The [eight-method table](qfo_fas_stratum_weights_20261001.md) and
[machine-readable report](qfo_fas_stratum_weights_20261001.json) separate
rounding from omission-induced changes to the native sample mixture. All
native benchmark scores, endpoints, defaults and earlier limitations remain
unchanged. This is not a new benchmark or a confidence interval.

## Sampling Contract

The retained scorer contains two strata: P eligible pairs with precomputed
scores and M annotated eligible pairs without those scores. It requests
c=min(M,9000) newly calculated pairs and k precomputed pairs, where k is the
rounded native floating-point calculation c*(P/(P+M))/(1-P/(P+M)), limited
by the available precomputed count. It shuffles the precomputed stratum and,
when M>9000, the newly calculated stratum, without setting a seed. The saved
score is the arithmetic mean of returned values, not a post-stratified mean.
All eight retained methods have M>9000 and request 9,000 new scores.

Write A and B for the saved precomputed and new-score means and r for the
number of returned new scores. All k requested precomputed scores return.
The three precomputed weights are:

```text
w = P/(P+M)          eligible-population composition
a = k/(k+c)          intended sample composition
b = k/(k+r)          saved sample composition
native mean = b*A + (1-b)*B
diagnostic = w*A + (1-w)*B
native minus diagnostic = (a-w)*(A-B) + (b-a)*(A-B)
```

The first component is rounding; the second is the mixture-weight effect
of score attrition, holding A and B fixed. The diagnostic is not an estimate
of the unknown eligible-population mean. In particular, omitted scores may
also change B; neither this decomposition nor its absolute weight-only
bound accounts for that change or for unsampled pairs.

Under hypothetical uniform sampling within both strata and complete score
return, the native mean's expectation has weights a and 1-a, not exactly
w and 1-w. If the rounded count is within half a unit of c*P/M, then
abs(a-w)<=M/(2*(P+M)*(k+c)). Exact rational arithmetic confirms that
half-unit condition in all eight retained rows. This does not validate
historical RNG state, pair inclusion probabilities or ignorable omissions.
Protein/family dependence prevents interpreting the pair-based native SEM
as paired method uncertainty. No new uncertainty estimate is admitted.

## Observed Arithmetic

OrthoMCL's intended precomputed weight is 91.300061%, but its 1,252 omitted
new scores raise the saved weight to 92.418564%. Holding the observed
stratum means fixed, the attrition component is +0.001318979 FAS units;
the rounding component is -0.000000033. The total is +0.001318947.
For full OrthoFinder, the corresponding total is +0.000085894; for
high-sensitivity and phylogenetic OrthoHMM it is +0.000008797 and
+0.000030716. These are numerical differences from a diagnostic, not
estimated biases relative to a known population truth. They do not establish
or overturn a full-population method advantage.

The previously retained intended-sample completion bounds and eligible-
population outer bounds address different missing-data questions. They are
not superseded. The actual omitted pair identities, representativeness and
dependence-aware FAS comparison uncertainty remain unresolved.

## Evidence And Validation

The new standard-library audit checks the fixed original order of all eight
methods, exact report pins, native request arithmetic, matching retained
population recounts and flags, full precomputed return, attrition totals and
saved endpoint reconstruction. It reads only two small retained reports;
no large lookup scan, database join, annotation, inference or FAS rescore is
repeated. The report records and rechecks both inputs and its source/helper.

All 115 focused sample/attrition/population/new-diagnostic tests pass. A
separate isolated standard-library Decimal calculation, without importing
the producer, checked all eight rows and four current file pins. Its largest
component discrepancy was 1.2271e-16; all native means remained identical.
This is independent aggregate arithmetic, not independent FAS scoring.

Reproduce into fresh output paths:

```sh
/usr/bin/python3 -B -m benchmark_tools.audit_fas_stratum_weights \
  --output /tmp/qfo_fas_stratum_weights.json \
  --table /tmp/qfo_fas_stratum_weights.md
```

The existing native container probe records `/benchmark/fas_benchmark.py`
with SHA256 `1045c57f4d0f4787bec3d1f0690799df63c68dcc67ccfd5a925c338e3d33661d`,
matching the frozen checkout source. That existing evidence is not a repair
of the separate original population-audit parser/stability provenance gap.
No historical provenance or controlled-timing status is upgraded here.
