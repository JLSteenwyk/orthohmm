# Retained Conditional Native FAS Result

## Actual Application

The [tested adapter](../report_native_fas_conditional.py) was pushed as c7200cc4
BEFORE applying the [frozen reporting plan](NATIVE_FAS_CONDITIONAL_REPORTING_PROTOCOL_20261010.md).
Focused/adjacent tests EXIT0: 151 passed in 2.50s. Actual retained-data
application EXIT0: all eight conditional ranges and all 28 contrasts computed.
No native benchmark score, annotation, RNG sample, inference or timing rerun.

| Method | Observed Native Z | Conditional 95% Range For Expected Theta |
| --- | --- | --- |
| orthohmm_high_sensitivity | 0.776585042 | [0.771332967, 0.780286193] |
| orthohmm_phylogeny_satellite_v2 | 0.762993312 | [0.761929649, 0.763906702] |
| orthofinder_3_1_5_full | 0.691422133 | [0.688429518, 0.694915139] |
| orthofinder_3_1_5_sequence_only | 0.561753084 | [0.546499945, 0.577724179] |
| sonicparanoid_2_0_9 | 0.736680187 | [0.732115268, 0.738512044] |
| proteinortho_6_3_6 | 0.813594864 | [0.813565174, 0.813689528] |
| fastoma_0_3_5 | 0.654366336 | [0.643583534, 0.662539506] |
| orthomcl_1_4 | 0.733752216 | [0.731706501, 0.735008398] |

Observed Z is preserved exactly in the [full JSON](native_fas_conditional_report_20261010_v1/report.json).
Decimals above are rounded. Theta is the expectation of the repeated native
post-attrition ratio under the stated fixed-population design, not the fixed
observed score. Unknown G/mu_G preclude an exact theta point value. Ranges
must not be plotted as biological error bars centered on observed Z.

The [complete generated table](native_fas_conditional_report_20261010_v1/summary.md)
retains all 28 left-minus-right differences, including negative ranges and
zero-overlapping comparisons. The single joint rectangle allocates alpha=0.05
across the count and mean components for eight methods; no additional cross-
method independence or favorable-comparison selection is imposed.

## Comparisons And Limits

On this conditional expected-FAS target only, high-sensitivity OrthoHMM minus
full OrthoFinder has range [7.641783, 9.185667] percentage points;
phylogenetic OrthoHMM minus full OrthoFinder has range
[6.701451, 7.547718] percentage points.
High sensitivity minus Proteinortho is unfavorable:
[-4.235656, -3.327898] percentage points.
Phylogenetic OrthoHMM is also below Proteinortho, and below high sensitivity.
The SonicParanoid/OrthoMCL conditional difference overlaps zero. These results
do not establish overall orthology accuracy, biological causality, superiority
on other metrics or generalization to new families/clades. Earlier negative F1
comparisons remain unchanged; no scoring endpoint or default was tuned here.

## Independent Check

The [actual command/terminal record](native_fas_conditional_execution_20261010_v1.json)
retains execution, tests and an independent 80-digit Decimal recurrence
readback. Before checking native values, that reader tests its relative-PMF
recurrence against all 15 exact rational validation targets (maximum error
1e-80). It then checks eight native integer count ranges, 32 endpoint/neighbor-
condition checks including support-boundary shortcuts, conditional-mean bounds,
E[k/(k+R)] weights, all 28 projections, preserved observed values, table values,
all direct source/input bindings and non-admission flags. It imports neither
the producer/kernels nor SciPy and normalizes the entire hypergeometric support.
Maximum floating/reference difference: 2.825517597671e-14;
readback EXIT0. This numerical check does not establish missing historical
source identity or the design assumptions themselves.

Report: 21615 bytes, SHA256
`83548e93afd5ffe0eeeef13250684dd65485b9efba1934d3377787659f7b3d7f`.
Generated table: 4709 bytes, SHA256
`55dc72d5f0da6c5432fcae6e8252919494621f5a4e51c4612027750ae59275ff`.
Only the new postprocessing driver/readback used one own CPU, four-GiB address-
space limit and one numerical thread; unrelated work was not altered. These
are resource limits, not measured peak-memory figures or isolated tool timings.

## Admission Scope

The analysis conditions on retained eligible populations/precomputed means,
uniform independent within-method stratum selections and fixed numeric-return
status/value on normal successful execution. No missing-at-random or biological
pair/family independence is assumed. The [native finite context check](NATIVE_FAS_CONTEXT_RESULT_20261010.md)
and [finite design validation](NATIVE_FAS_SAMPLING_VALIDATION_RESULT_20261010.md)
support this explicitly conditional analysis, not every historical annotation,
file race, batch crash, worker assignment or ideal historical RNG state.

Historical parser and database hash bindings remain unestablished/false.
Missing historical pair identities are not reconstructed. The retained-data
conditional ranges are computed; unconditional historical interval admission,
biological generalization, other-endpoint uncertainty and full publication
readiness remain false. Old reports/flags and all observed benchmark scores
are preserved. The full publication goal remains incomplete.
