# Retained QfO FAS Sample Audit

The saved pair scores reproduce all eight historical comparator means and
standard errors and all four recovered-stage means and standard errors within
1e-12. No historical score or method setting was changed. Machine-readable
evidence, input hashes, sample protein degrees and all 66 pairwise sample
overlap counts are in `qfo_fas_sample_audit_20260917.json`.

[Generated table](QFO_FAS_SAMPLE_TABLE_20260917.md) reports all twelve samples.
Coverage spans 0.0067% for the OrthoFinder MCL checkpoint to 31.3477% for
ProteinOrtho. FastOMA used a supplied tree. Recovered stages are separate
from the historical comparator runs and are not phylogenetic ablations.

## Native Statistic and Sampling

The inspected retained `fas_benchmark.py` averages precomputed directional
values into one score per canonical accession pair. It counts precomputed
pairs plus missing pairs with both proteins annotated as `NR_ORTHOLOGS`;
this is neither the scored sample size nor all predicted relations. It drops
accession aliases containing underscores and reports unannotated pairs
separately in logs.

When missing scores exist, it samples the precomputed stratum to maintain its
original fraction relative to at most 9,000 newly computed pairs. Both
shuffles use Python's random module without a seed set in this scorer.
Available newly computed results are included; missing results are skipped.
Failure or incomplete new scoring can therefore change the realized stratum
mixture. Historical RNG state and actual stratum-specific inclusion rates
are not reconstructed by this audit. The cap applies to new calculations,
not the total sample, explaining samples above 9,000 pairs.

The reported mean is the arithmetic mean of saved scores. Native stderr is
sample standard deviation (`ddof=1`) divided by square root of sample size.
Every retained sample contains proteins used in multiple pairs. Native SEMs
must not be interpreted as family-aware uncertainty for method differences.
Small sampling fractions alone do not establish bias, but realized coverage,
missing-score handling, dependence and seed provenance need to be reported.

## Scope and Next Steps

This is an arithmetic and sample-structure audit, not validation of the
underlying FAS implementation, annotations, prediction-to-pair conversion or
sampling representativeness. `reported_eligible_pairs` is checked for type
and consistency with sample size, not independently recounted from databases.
Pair overlap does not define independent resampling units. A defensible
comparison interval still requires a dependence-aware analysis and explicit
treatment of method-specific sampling; bootstrapping saved pairs as IID is
not an adequate substitute. Retain native scores, do not replace them with a
shared-pair-only endpoint selected after seeing results.

Reproduce from retained raw inputs:

```bash
python benchmark_tools/audit_qfo_fas_samples.py --output NEW_AUDIT.json --table NEW_TABLE.md
python -m pytest -q tests/unit/test_audit_qfo_fas_samples.py
```
