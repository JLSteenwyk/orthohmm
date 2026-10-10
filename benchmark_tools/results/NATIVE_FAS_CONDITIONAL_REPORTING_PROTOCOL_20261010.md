# Conditional Native FAS Reporting Protocol

Prospective original4.1 follow-through. No historical eight-method interval has
been computed or viewed. This plan uses the actual native sampling statistic
and retained data, not a new complete-case or post-stratified endpoint.

## Target And Claim Boundary

Report confidence sets for the expected repeated post-attrition native FAS
ratio under the fixed-population design in
[the derivation](NATIVE_FAS_SAMPLING_DESIGN_PROTOCOL_20261010.md).
The [finite coverage validation](NATIVE_FAS_SAMPLING_VALIDATION_RESULT_20261010.md)
and [controlled native contexts](NATIVE_FAS_CONTEXT_RESULT_20261010.md) support
this conditional route, not unconditional historical certification.

Require the explicitly stated assumptions: the retained eligible population
and precomputed mean describe the conditioned population; both strata are
sampled uniformly without replacement and independently within a method;
successful native scoring gives fixed per-pair numeric return status/value
independent of sampled companions/order/worker assignment. Missingness may
depend on a pair's own features and score. No missing-at-random assumption,
biological pair/family independence or cross-method independence is made.

The native observed mean Z and expected statistic theta are distinct. Preserve
every observed mean exactly in a separate descriptive column. Unknown G and
mu_G prevent an exact point value for theta: report its confidence range, not
an expected-count substitution, midpoint estimate or relabeled observed mean.
Do not attach these ranges as biological error bars to observed Z. Do not
describe them as new-family/clade uncertainty, overall accuracy superiority,
an interval for another QfO endpoint or independent biological confirmation.

## Frozen Inputs

| Input | SHA256 |
| --- | --- |
| qfo_fas_sample_attrition_20260928.json | 1f6ddfd5aadf9f5be8147164b23904d2e5db1ffee69a6fc18b654a54f3f25d36 |
| qfo_fas_population_completed_22383/report.json | 05f1ab3f7f0ffed8ddcffe2a6bfd4ce26ef46b974dbb946cb60883857803eeaf |
| native_fas_sampling_validation_20261010_v1.json | 236ffb82887d816904d08df0e5b350a358badd8e9f5f6fc088004cd0d021ac57 |
| native_fas_context_probe_20261010_v2.json | 845b4baea528f7633830e5a26ec58267da95165d8c4f79f2a0112c6063826f33 |
| ../native_fas_sampling_interval.py | 803ff0f455213e6a55b83e30c586ab1c57566830f8f877c916b78e84bf3658f9 |
| ../audit_fas_stratum_weights.py | 8712b7e36fc756b61c161847493f7c65b7b376d5ef675c6c803f52ab21d5fbc3 |

Use existing `audit_fas_stratum_weights.panel` for unchanged identity/count/
sample checks, rather than copying its validation or rerunning its producer.
Keep the eight methods in its frozen METHODS order. P/M and mu_P come from
the population report. k/c/r and the returned-score sample mean come directly
from the attrition report's original strata. Do not substitute the sampled
precomputed mean for the full precomputed population mean. Preserve omitted
return counts and the original saved mean. No large database/lookup scan,
raw rescore, new annotation or new historical RNG sample is needed.

## Numerical Analysis

Apply the tested frozen kernel with alpha=0.05 and component error=0.05/16.
Invert count tails, retain the full conditional-mean uncertainty and project
using E[k/(k+R)]. Report G bounds and the expected-ratio bounds for every method.
List all 28 contrasts in the fixed method order, sign left-minus-right,
including intervals overlapping zero and differences unfavorable to OrthoHMM.
The shared rectangle supplies simultaneous coverage; do not add another
independence assumption or select only comparisons that exclude zero.

Validate numerical support, finite values, bounds ordering and source/input
identities. Do not weaken the kernel's tail or mass checks to force output.
If a method fails, preserve that new postprocessing failure and leave its
dependent contrasts unavailable. Do not relaunch historical scoring or fill
missing values. Any supported numerical correction must be separately named,
tested and frozen before use; retain the original attempt unchanged.

## Provenance And Admission

The population audit reports historical parser hash identity and historical
database hash binding as unestablished. Those remain false and prominent.
Agreement of retained scalar summaries does not establish them. Historical
RNG states and omitted-pair identities are unavailable; finite context checks
do not certify every historical annotation, file race, crash or worker choice.

A successful report may state conditional design intervals were computed
under the assumptions above. It may not assert unconditional historical
interval admission or biological generalization, change old admission flags,
replace published point scores or declare the full publication goal complete.
Preserve old evidence and its scopes. Keep other endpoint uncertainty, failed
factorial cells, missing TreeFam originals and full-goal gaps explicit.

Implement a thin reporting adapter with focused tests; commit/push source
before its first prospective application to these retained values. Use fresh
outputs, observe the same handle to terminal and integrate actual results,
including failures, without replaying a completed producer.
