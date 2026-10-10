# Native FAS Finite-Sampling Validation Result

## Actual Outcome

The [prospectively frozen driver](../validate_native_fas_sampling.py) completed
EXIT0 after source/tests were pushed as e5fd810d. The
[full result](native_fas_sampling_validation_20261010_v1.json) contains all
15 predeclared cells. This validates the conditional numerical construction
in the [native design protocol](NATIVE_FAS_SAMPLING_DESIGN_PROTOCOL_20261010.md),
not historical benchmark intervals or biological generalization.

The expanded focused/adjacent tests passed 132 tests in 2.30s before execution.
The [actual execution and readback record](native_fas_sampling_validation_execution_20261010_v1.json)
retains commands, yielded observations and terminal results. No native scorer,
new historical sample, annotation, sequence inference, timing panel, RNG draw
or omitted-probability-tail approximation was used.

| Validation Cell | M / c / G | Expected Native Ratio | Target Miss Probability |
| --- | --- | --- | --- |
| zero_returns | 8 / 4 / 0 | 0.5000000000 | 0 |
| singleton_return | 1 / 1 / 1 | 0.5500000000 | 0 |
| zero_draws | 8 / 0 / 4 | 0.5000000000 | 0 |
| complete_draws | 8 / 8 / 5 | 0.5227272727 | 0 |
| random_denominator_example | 4 / 2 / 2 | 0.4458333333 | 0 |
| balanced_partial | 128 / 64 / 64 | 0.5000000000 | 3.272396e-7 |
| high_pre_low_return | 128 / 64 / 80 | 0.3203661540 | 1.865357e-12 |
| low_pre_high_return | 128 / 64 / 80 | 0.6796338460 | 1.865357e-12 |
| rare_partial_returns | 128 / 32 / 2 | 0.5000000000 | 0 |
| all_return_partial | 128 / 64 / 128 | 0.5000000000 | 1.811388e-7 |
| native_cap_sparse | 12000 / 9000 / 2 | 0.5159375781 | 0 |
| native_cap_dense | 12000 / 9000 / 11999 | 0.2501249688 | 0 |
| native_cap_no_returns | 12000 / 9000 / 0 | 0.5000000000 | 0 |
| native_cap_all_returns | 12000 / 9000 / 12000 | 0.7499166944 | 0 |
| unequal_strata | 5 / 3 / 2 | 0.7424038462 | 0 |

M is the missing-lookup population, c the requested sample and G the fixed
numeric-return population. Four native_cap cells use the literal c=9000 with
M=12000 and positive native k=3. Other cells are explicit small-population
laboratory caps or boundary cases, not historical benchmark executions.
Nonzero miss probabilities are shown even when rounded coverage would look
like 1. Full rational probabilities and targets remain in the JSON.

## Target And Coverage Checks

Uniform finite subsets are collapsed into category counts with exact integer
multiplicities. Both strata's subset scores are independently averaged through
the actual random-denominator ratio. For the diagnostic example the exact
answer is 107/240, not the expected-count substitution 9/20. The literal-cap
sparse case likewise has a nonzero substitution error. These are examples of
why the random denominator must remain in the target, not new endpoints.

Count endpoints are independently checked by rational hypergeometric tails,
including exclusion of adjacent integer values. Mean bounds and rectangle
corners use an 80-digit Decimal reference with a stated 1e-60 endpoint-equality
margin. Coverage is the exact rational mass of reference-classified outcomes;
floating-kernel agreement is checked separately with 1e-12 tolerance. This is
not a claim of arbitrary-precision output from the floating kernel.

Minimum count coverage: 0.997482469652.
Minimum conditional-mean and target coverage: 0.999999672760.
Maximum floating/reference error, including differences:
3.330669073875e-16.
The component error is 1/320; each method has the union-bound error limit 1/160.

The eight selected validation cells share a uniform rank coupling, with
alternate outcome order reversed. Independent uniform precomputed selections
are retained within each method. Across 2695 exact coupling segments, joint
mean coverage is 0.999999491620 and joint coverage of
all 28 differences is 0.999999934563. Swapping
methods negates and reverses the endpoints. The summed component failure bound
is 0.006053906575, below alpha=0.05.
The union-bound argument permits arbitrary cross-method dependence; the
specific coupling checks implementation but does not enumerate every coupling.
These cells are not eight real tools or independent biological replicates.

## Rejected Control And Remaining Scope

An explicit violation has A return only when B accompanies it. Its actual
expected ratio is 22/45; pretending singleton returns define a fixed population
gives 8/15. The fixed-return assumption is rejected. No coverage or native
admission is claimed for that control.

Separate readback recomputes all 15 rational targets from hypergeometric
count probabilities and conditional native ratios, checks the plugin
contrasts, reported component/joint bounds, control, all source/protocol
identities and false non-admission flags. Readback EXIT0. It is not a second
independent reconstruction of every coverage classification. Result identity:
17766 bytes, SHA256
`236ffb82887d816904d08df0e5b350a358badd8e9f5f6fc088004cd0d021ac57`.

The [native controlled context check](NATIVE_FAS_CONTEXT_RESULT_20261010.md)
supports fixed pair outcomes only for its successful finite fixtures. This
validation assumes uniform within-stratum sampling, independent strata within
a method, fixed return status/value and a known precomputed population mean.
It allows score-dependent omissions and assumes no biological gene/family
independence or missing-at-random mechanism. It does not repair historical
parser/database hash-binding gaps, establish all historical worker/file
behavior, recover missing pair identities or prove ideal historical RNG states.

No historical eight-method intervals have been viewed or admitted. Their
reporting plan and target/provenance qualifications must be frozen first.
The expected repeated native ratio is distinct from the observed fixed mean,
full eligible-population accuracy and new-family/clade generalization. Other
native endpoint uncertainty and full publication readiness remain unresolved.
