# Native Profile-Refinement Uncertainty

This chronological supplement adds P1C0R1 to the
[four-cell native accuracy table](native_qfo_scientific_scores_20261007_v3/scores.md)
and the retained paired-family analysis. Earlier three-cell manuscript,
figures and archives remain historical and unchanged. This is conditional
downstream profile refinement at C0/R1; initial HMM search is on in every cell.
It does not isolate the value of all HMMs or compare directly with OrthoFinder.

## New Evidence

Before selected execution, adapters/tests/protocol were committed and pushed
at a9545d27. The new
[native10 count audit](native10_allocated_swiss_counts_20261007_v1.json), SHA256
`d5bafed830bf93e46b810728df14ffa34185d58ffe3eac78a693b0159d3d0a40`,
reads only the newly independently admitted P1C0R1 raw output, using the
unchanged original count/native-endpoint kernels. All18 family records,
members, reference truth universe and aggregate match the retained corrected
counts exactly. Native decimal SwissTrees F1 is0.7863721859036304;
count-derived F1 is0.7863721816702856. This tiny serialization difference
does not define a new endpoint. Previous raw audits were not recounted.

The [four-cell interval binding](native_qfo_four_cell_swiss_uncertainty_20261007_v1.json),
SHA256 `e8f554a178841da24a1892a7e42f8680f9fe8e4827dec6d6bfe4643cdb890bac`,
reuses the original frozen projector and100,000 paired18-family draws,
seed20260922, with all42 endpoints retained in adjustment. Four cells bind;
three conditional contrasts now match (P_at_C0_R1, C_at_P0_R0, R_at_P0_C0),
and11 remain unavailable with null metrics. No interaction is supported yet.
Only complete family-record/aggregate equality permits interval reuse.

## Independent Readback

Reader source was committed/pushed at f12e72fb before one selected execution.
The [readback](native_qfo_profile_swiss_readback_20261007_v1.json), SHA256
`3f62b487014cb9ba99ea36b5adc15535f62c7e39c9a99f3cb4eb8c537523005d`,
uses the separate existing CSV/gzip parser and exact Fraction arithmetic,
not primary count/bootstrap functions. It checks36 native family records,
32,295 raw rows and10,765 P1C0R1 scored pair labels against retained evidence,
with identical truth/member universes. It reconstructs macro precision and
recall first, their harmonic mean second, conditional differences and family
signs. Both prior matched contrasts are exactly unchanged. Interval numbers
are checked against the frozen bootstrap, not independently recalculated or
certified. It does not prove identical whole partitions or other endpoints.

Count audit/binder/readback each executed once and exited0, in3.95/3.47/2.17s
with517,260/517,264/26,944KiB maximum RSS and zero swaps. Exact commands/outcomes
are in [audit execution](native10_swiss_audit_execution_20261007_v1.json),
[binding execution](native10_swiss_bind_execution_20261007_v1.json) and
[readback execution](native10_swiss_readback_execution_20261007_v1.json).
These are read-only analysis costs, not inference or isolated tool efficiency.
All selected source/report bytes are now bound; do not edit or rerun them
merely because another continuation starts.

## Conditional Result

| Metric | Difference (pp) | Nominal 95% Paired CI | 42-Endpoint Adjusted CI | Family Wins/Ties/Losses |
| --- | ---: | --- | --- | --- |
| SwissTrees F1 | -0.320215 | [-1.256479, 0.262365] | [-2.112886, 0.513227] | 1/16/1 |
| SwissTrees precision | -0.010056 | [-0.793651, 0.763483] | [-1.322751, 1.272471] | 1/16/1 |
| SwissTrees recall | -0.462963 | [-1.388889, 0.000000] | [-2.314815, 0.000000] | 0/17/1 |

Profile refinement has a slightly lower F1 point estimate with an adjusted
interval including zero. This is not equivalence, significant improvement,
general inferiority or default-tuning evidence. Sixteen families have unchanged
per-family metrics. CASP improves in F1/precision without changed recall;
GH14 worsens in all three. These signs are observed, not an established
mechanistic explanation; family tracing remains a separate task.

The analysis is development-exposed and conditional on18 families, approximate
percentile coverage and exchangeability assumptions. Other QfO endpoints and
the project-defined secondary mean do not inherit these intervals. R-on
resolved-pair versus R-off group-clique semantics remain explicit. The
recovered reference's native timing remains failed, null and ineligible; no
timing repair or scientific timing admission follows. The factorial and
publication package remain incomplete, with failure9 retained and genuinely
unrun score identities11/12 still unavailable at this result checkpoint.

## Validation And Reproduction

Expanded count/binding suite passes193 tests in8.99s; independent-reader joined
suite passes196 in9.93s, no errors/failures/skips. The earlier two failed fixture
checks remain retained; only their overwritten synthetic snapshot path was
fixed, without weakening production provenance. Tests use invented handoff
metadata, not biological truth or a claim of production runtime closure.
[Count/binding protocol](ALLOCATED_NATIVE_SWISS_UNCERTAINTY_PROTOCOL_20261007.md)
and [readback protocol](ALLOCATED_NATIVE_PROFILE_READBACK_PROTOCOL_20261007.md)
define execution scope; the exact retained commands identify actual outputs.
For genuine reproduction use fresh destinations and preserve existing files.
No new native inference, conversion, scoring, FAS sample or bootstrap is needed.

Timing measurements were collected on a shared Threadripper while other
analyses were running. Competition for CPU, memory bandwidth and I/O may
have affected elapsed times, with an unknown and potentially tool-dependent
impact. These are observed shared-host timings, not estimates of isolated
performance.
