# Native SwissTrees Domain-Bin Results

[Fixed protocol](NATIVE_QFO_SWISS_DOMAIN_STRATA_PROTOCOL_20261006.md) pushed at
`642ed16d` before this native projection. Exporter, independent reader and
33 new focused tests pushed at `880acc47` before actual computation. No method
tuning, new inference, new scoring admission or historical interval transfer.
The [score table](native_qfo_swiss_domain_strata_20261006_v1/TABLE.md),
[structured report](native_qfo_swiss_domain_strata_20261006_v1/report.json) and
[independent readback](native_qfo_swiss_domain_strata_readback_20261006_v1.json)
are new native evidence, not relabeled selected-default results.

## Point Estimates

F1 percentages; changes are percentage points for `p0_c0_r1` minus `p0_c0_r0`.
Initial HMM search remains on; profiles and candidate expansion are off in
both. R0 uses group-clique predictions, R1 resolved native pairs. All rows
use macro family precision/recall followed by harmonic F1, with the original
TP/FP/FN divided-by-two plus-one prior, not pooled pairs or mean family F1.

| Fixed Bin | Families | R0 F1 | R1 F1 | F1 Change | Precision Change | Recall Change |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| All | 18 | 68.918 | 78.957 | +10.039 | +30.521 | -6.533 |
| Median Pfam types < 2 | 12 | 71.006 | 75.897 | +4.891 | +24.288 | -9.169 |
| Median Pfam types >= 2 | 6 | 63.986 | 84.932 | +20.946 | +42.987 | -1.262 |
| Repeated-type fraction < 1/4 | 15 | 69.545 | 79.749 | +10.204 | +29.147 | -5.181 |
| Repeated-type fraction >= 1/4 | 3 | 65.595 | 74.731 | +9.136 | +37.390 | -13.296 |

The larger F1 point gain in the six higher-type families accompanies a larger
precision gain and smaller recall loss. The three-family repeated-type bin
instead loses 13.296 recall points despite higher F1. Neither observation
establishes a domain-driven cause or statistically demonstrated interaction.
Both splits are unchanged September 17 definitions and overlap; all 18
families and all 36 cell/family count records remain in the report.

## Verification And Scope

Both actual executions exit zero. Each reads the original native raw files
and all 49 original annotation JSONs containing selected accessions, checking
433,906,275 annotation bytes against original inventory hashes. Independent
annotation parsing and alternative macro/harmonic arithmetic agree on all
563 exact accession features, including five zero-hit annotations, 21,530
raw relations, ten score rows and five differences. The remaining 29 annotation
source identities are inherited only, not newly checked. Direct checks do
not constitute a new transitive admission or prove uninterrupted integrity.

The final joined [JUnit receipt](native_qfo_swiss_domain_strata_20261006_v1/native_qfo_swiss_domain_strata_joined_tests_20261006_v1.xml)
reports 92 tests passed in 1.68s, no failures/errors/skips: 33 new cases and
existing sequence-stratum, annotation extraction and historical-bin contracts.
The earlier 33-case receipt is retained. Synthetic refusals cover missing,
extra/shared accessions, changed summaries/bins, altered or incomplete counts,
wrong statistics/semantics, duplicate rows/instances, invalid coordinates and
existing/dangling output paths. Actual export/readback success is separately
recorded; fixture success is not scientific replication.

Original Python 3.10.13 invocation uses the retained review environment, with
Python/site/preload/library overrides removed and bytecode disabled; the
independent reader also uses `-I -S`. No scientific environment changes.
[Export receipt](native_qfo_swiss_domain_strata_20261006_v1/native_qfo_swiss_domain_strata_export_20261006_v1.time.txt)
reports 32.71s/396,028KiB maximum process RSS;
[reader receipt](native_qfo_swiss_domain_strata_20261006_v1/native_qfo_swiss_domain_strata_readback_20261006_v1.time.txt)
32.10s/372,148KiB. Both exit zero and report zero swaps. These are shared-host
postprocessing measurements, not inference costs. CPU, memory-bandwidth and
I/O contention have unknown, potentially tool-dependent effects; no isolated
speed ranking. At launch 671,965,057,024 available RAM bytes and only
25,817,088 free swap bytes were observed; unrelated workloads remain untouched.

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| `export_native_qfo_swiss_domain_strata.py` | 12525 | `ccc979a24eeb78965081e72a78da9fa9fc0eaca50ead56c59e3ab0c8e5230ae0` |
| `readback_native_qfo_swiss_domain_strata.py` | 12375 | `a6048b7518c26d7c86d99cd049fc8d625b04be3dbea48311875f1baf22788e66` |
| `native_qfo_swiss_domain_strata_20261006_v1/report.json` | 80268 | `d3235c712e711f005b3182a3621c0079f8a85c3c7abd5db016d5e8ea95884cd6` |
| `native_qfo_swiss_domain_strata_readback_20261006_v1.json` | 26139 | `9e77b67ad97e6824e9442691be7a05ba1552d62ffe696ced6034a996ede1945d` |

## Remaining Requirements

No new confidence intervals or bootstrap draws: historical domain intervals
belong to different selected-default predictions. Small curated, development-
exposed bins are not independent confirmation. Pfam counts are not complete
validated domain architectures, fragment/domain-loss truth, calibrated
divergence or ancestral duplication histories. Shared FAS annotations are not
independent FAS validation. This advances native error-stratum evidence without
closing wider goal 4.3 or other-endpoint uncertainty.

R1's failed original timing remains ineligible; native scores stay admitted
separately. No original job restarted or unfinished result inspected. At last
check original 22444 RUNNING 10:23:49; 22445/22450/22451/22452 dependency-pending.
Await their original terminal gates before new admission or the next native
identity. Five native score cells remain unavailable. Historical supplement,
main text, archive and frozen execution/scientific protocol bytes unchanged;
the full publication goal stays active and completion remains unproven.
