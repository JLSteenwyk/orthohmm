# Corrected FAS Sample Audit

The [current eight-method table](QFO_CORRECTED_FAS_SAMPLE_TABLE_20260928.md)
reproduces corrected-release FAS means and native standard errors from saved
pair scores, within absolute tolerance 1e-12. Unlike the historical audit,
these samples belong to the corrected v7 comparison, including recovered
OrthoMCL. Each raw file is bound to its assessment execution-output checksum;
the aggregate, participant identity and eligible-pair count match the admitted
comparison. All checked records are rehashed after reading.

Scored fractions range from 0.0067% for the OrthoFinder MCL checkpoint to
58.4605% for ProteinOrtho. High-sensitivity and phylogenetic OrthoHMM score
38,158 and 175,357 pairs (0.4224% and 2.9424% of reported eligible pairs);
full OrthoFinder scores 53,367 (0.3754%). All eight samples reuse proteins
across pairs. These observations describe the samples, not evidence that
sampling fractions alone cause bias or explain a method's score advantage.

The [native sampling review](QFO_FAS_SAMPLE_AUDIT_20260917.md) remains relevant:
NR_ORTHOLOGS is not the scored sample size, and the newly calculated-pair cap
does not cap total saved pairs. This audit does not reconstruct random states,
stratum inclusion probabilities or missing-score mechanisms. It does not
independently recount eligible predictions or validate the underlying FAS
algorithm and annotations. Reproducing native pair-IID SEMs does not validate
them as family-aware uncertainty for a method difference. No confidence
interval, endpoint change or new accuracy claim is admitted.

The [machine-readable audit](qfo_corrected_fas_sample_audit_20260928.json)
records sample protein degrees, source hashes and all checked inputs.
Reproduce into fresh paths:

```sh
python -B -m benchmark_tools.audit_corrected_fas_samples \
  --output NEW_AUDIT.json --table NEW_TABLE.md
python -B -m pytest -q tests/unit/test_audit_corrected_fas_samples.py \
  tests/unit/test_audit_qfo_fas_samples.py
```
