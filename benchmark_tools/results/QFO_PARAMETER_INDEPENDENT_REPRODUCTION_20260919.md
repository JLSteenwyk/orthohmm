# Independent Parameter Uncertainty Reproduction

`benchmark_tools/reproduce_qfo_parameter_uncertainty.py` independently
reconstructs the frozen six-contrast parameter analysis from admitted raw
SwissTrees family counts. It does not import the production bootstrap,
score calculator or count auditor. It computes the native smoothed family
precision/recall from TP/FP/FN, then recomputes macro precision/recall and
their harmonic mean in each paired family resample. Weighted family sums
replace the production matrix multiplication.

The controls remain 100,000 draws, PCG64 seed 20260925, 18 families, six
candidate-minus-control contrasts and 18-endpoint Bonferroni adjustment.
Missing arms do not reduce multiplicity or become zeros. The verifier
checks points, nominal and adjusted intervals, family differences and
wins/ties/losses within absolute tolerance 1e-12; it also checks availability,
reference membership/truth totals and completeness flags. File records
are rehashed before and after numerical reconstruction. Existing outputs
are not overwritten.

## Verification

63 focused tests passed in 11.07 seconds across the independent checker,
production kernel and provenance-bound runner suites. Tests cover all 18
synthetic endpoints, partial panels, unavailable baseline, control-only
panels, corrupted counts/statistics/intervals, array-shape mismatch, NaN,
changed file records and output overwrite prevention.

The actual control-only result from job 21986 was checked with:

```sh
python -B -m benchmark_tools.reproduce_qfo_parameter_uncertainty \
  --results benchmark_tools/results/qfo_parameter_uncertainty_control_21986.json \
  --results-sha256 c930b10fde59c799bc8d7f297a958320d07e5b45c934433b6afe49cafc837c47 \
  --output benchmark_tools/results/qfo_parameter_uncertainty_control_reproduction_21986.json
```

The [retained report](qfo_parameter_uncertainty_control_reproduction_21986.json)
records successful numerical reproduction under NumPy 2.2.6 and source
SHA-256 `6b337c936bf1cf6896fea06de6e1924092c1dc6e5b3871ad621edbdc2138ddf8`.
There are **zero estimated endpoints**: only the real control is admitted.
Its point estimates and family statistics agree, and all six unavailable
contrasts retain null results. This is integration evidence, not observed
parameter stability or uncertainty for variants.

## Limits And Next Step

This implementation is independent arithmetic, not an independent raw-data,
scoring or native-inference audit. Both implementations use NumPy PCG64 and
linear quantiles, so agreement cannot validate those shared primitives or
the statistical assumptions. The analysis remains development-exposed;
18 families may not be exchangeable. No defaults or endpoints were changed.

When genuine variant admissions arrive, freeze a new inventory, execute
the existing parameter uncertainty runner, and use this checker against
the resulting report and its exact checksum. Do not overwrite the
control-only integration evidence or infer equivalence from null intervals.
Parameter robustness and publication readiness remain incomplete.

## Replacement-Executor Bindings (2026-09-23)

All four replacement candidate-threshold inference runs, 22034_0..3, have
completed 0:0. Their 22035 native-validation array has started; conversion,
scoring and independent score admission remain queued. CPM replay 22081_0
is still running. No variant score is inferred from native completion.

The current uncertainty runner still pinned superseded executor identities;
five focused tests reproduced its fail-closed rejection. Reviewed Git diffs
show that the relevant scorer/admitter changes bind replacement job IDs,
executor revisions and v2/v3 paths. CPM context keeps control v1 and routes
replacement variants to v3. Updated only these source bindings in the count
auditor and uncertainty runner, preserving the frozen protocol, plan, counts,
statistics, bootstrap kernel, seed, multiplicity and historical executors.
A source-binding regression test now checks the two admitter identities too.

All 117 focused runner, count-audit, bootstrap and independent-reproduction
tests passed. A fresh local control-only integration run against the retained
raw evidence produced the [new report](qfo_parameter_uncertainty_rebound_control_20260923.json),
441,885 bytes, SHA-256
0a6e9b2a8001b207b94b45e25f228dddcf810d2168d97602cbe48b72c807ce8b.
Control point estimates and all comparison records equal the historical
21986 result. The [independent reproduction](qfo_parameter_uncertainty_rebound_reproduction_20260923.json)
also passes at 1e-12. There are still zero estimable variant endpoints,
complete_panel=false and uncertainty_admitted=false. This is repaired
integration evidence, not completed robustness or new scientific results.

The run used the same control-only inventory SHA-256
2cad3c5d168c49b4337668f14559f001dc1b7d7cd18b83f69668a208260e88c9,
baseline qfo_swiss_counts_20260917.json, frozen neighborhood plan/protocol,
and single-thread OMP/OpenBLAS/MKL settings. No scientific job was resubmitted.
When real admissions arrive, retain all six planned variants in a new frozen
inventory and execute the updated runner followed by numerical reproduction.
