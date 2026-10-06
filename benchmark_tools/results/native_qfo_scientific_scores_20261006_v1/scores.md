# Native QfO Scientific Results

| Index | Cell | Accuracy status | Measurement status | VGNC F1 | SwissTrees F1 | TreeFam-A F1 | GO similarity | EC similarity | FAS | Secondary mean | Submitted pairs | Input accessions | Relation accessions | Relation coverage | Prediction semantics |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 6 | p0_c0_r0 | supplied_native_admission | successful_native_terminal_review | 0.666834 | 0.689184 | 0.605404 | 0.472119 | 0.932114 | 0.774945 | 0.690100 | 9009082 | 984137 | 546450 | 0.555258 | cross-species group-derived clique pairs |
| 7 | p0_c0_r1 | supplied_recovered_scientific_admission | failed_timing_scientific_outputs_recovered | 0.898185 | 0.789574 | 0.602508 | 0.490260 | 0.967702 | 0.785007 | 0.755539 | 5113820 | 984137 | 542365 | 0.551107 | native phylogenetically inferred pairs |
| 8 | p0_c1_r0 | no_supplied_native_admission | no_supplied_accuracy_admission | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | cross-species group-derived clique pairs |
| 9 | p0_c1_r1 | no_supplied_native_admission | no_supplied_accuracy_admission | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | native phylogenetically inferred pairs |
| 10 | p1_c0_r1 | no_supplied_native_admission | no_supplied_accuracy_admission | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | native phylogenetically inferred pairs |
| 11 | p1_c1_r0 | no_supplied_native_admission | no_supplied_accuracy_admission | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | cross-species group-derived clique pairs |
| 12 | p1_c1_r1 | no_supplied_native_admission | no_supplied_accuracy_admission | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | Unavailable | native phylogenetically inferred pairs |

Timing measurements were collected on a shared Threadripper while other analyses were running. Competition for CPU, memory bandwidth and I/O may have affected elapsed times, with an unknown and potentially tool-dependent impact. These are observed shared-host timings, not estimates of isolated performance.

- Includes supplied successful-measurement and explicit recovered-science admissions, not live or cached scores.
- Recovered accuracy does not admit failed native timing; resource values remain null and eligibility false.
- The seven fresh identities exclude reused P1C0R0; historical/cached scores are not substituted.
- Initial HMM search remains on; R-off group cliques and R-on resolved pairs have different semantics.
- VGNC/SwissTrees/TreeFam-A are F1; GO/EC similarity and FAS are not F1.
- The six-metric mean is a project-defined secondary summary, not official QfO F1.
- Relation coverage uses all inference inputs, not reference-family coverage or accuracy.
- Direct report checks and endpoint arithmetic do not repeat transitive raw admission or biological validation.
- FAS sample/population limits and native SEM semantics are retained; no paired uncertainty or superiority claim.
- No inference, scoring, admission, job launch, timing correction or isolated-performance ranking.
- These are direct report/source/arithmetic checks, not a repeated transitive raw scientific or score admission.
