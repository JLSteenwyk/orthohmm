# FAS Requested-Sample Attrition

The [generated eight-method table](QFO_FAS_SAMPLE_ATTRITION_TABLE_20260928.md)
quantifies requested new FAS calculations without saved scores. This extends
the corrected sample arithmetic audit; it does not replace native endpoints
or supply the missing family-aware uncertainty analysis.

All eight native tasks requested 9,000 new scores. Missing counts are 7 for
high-sensitivity OrthoHMM, 27 for satellite_v2, 25 for full OrthoFinder, 36 for
its MCL checkpoint, 6 for SonicParanoid, 1 for ProteinOrtho, 49 for FastOMA and
1,252 for recovered OrthoMCL. The entire requested precomputed stratum is
present for every method. Counts sum to each admitted saved sample size and
the logged eligible strata sum to the admitted eligible count.

OrthoMCL's missing new scores represent 13.91% of its requested new stratum,
not 13.91% of its overall intended sample. They are distinct from sequence-
specific BLAST failures. Neither their identities nor their causes are
established here. Native code silently skips pairs absent from the final
score lookup; a successful task exit does not establish complete new scoring.

## Conditional Sensitivity Bounds

Let S be the sum of saved scores, n their count and d the requested new scores
without saved values. Assuming each omitted value would lie in [0,1], the
arithmetic mean of the intended sampled pairs lies in
`[S / (n + d), (S + d) / (n + d)]`. These are deterministic completion bounds,
not confidence intervals. They require neither independence nor a missing-at-
random assumption, but condition on the intended sample and the bounded-score
assumption. They do not cover unsampled eligible pairs or identify a corrected
population mean.

For OrthoMCL, this gives 0.724872 to 0.736975, compared with the saved mean
0.733752. The corresponding bounds are 0.776443 to 0.776626 for high-sensitivity
OrthoHMM, 0.762876 to 0.763030 for satellite_v2 and 0.691098 to 0.691567 for full
OrthoFinder. Completing only these omitted sampled scores cannot reverse the
two OrthoHMM-versus-full-OrthoFinder FAS sample-mean orderings. This narrow
observation does not address sampling error, representativeness, family
dependence, algorithm validity or general superiority.

## Binding And Validation

The [machine-readable result](qfo_fas_sample_attrition_20260928.json) retains
input hashes, stratum counts and full-precision bounds. The generated
table rounds to six decimal places, so ProteinOrtho's nonzero bound width is
not visible there. All 49 focused tests pass, and a second full audit
reproduces the complete JSON and generated table exactly. The prior corrected
sample audit's checked records are verified before and after this analysis.
Each native task is identified by its parsed `fas_benchmark.py` command and
participant, and its raw output must match the historically pinned published
raw file byte for byte. Downstream tasks share the results directory and are
explicitly excluded. Task logs and commands are newly hashed; this is not an
independent historical hash pin for their contents.

Native output lists precomputed scores before new scores. Splitting raw rows
by logged counts reproduces both six-decimal logged means, and the complete
sample reproduces the admitted mean within 1e-12. This does not independently
reclassify pairs using the full precomputed lookup. Eligible population counts
remain logged values, not independently recounted predictions. Historical
random states and omitted-pair identities remain unavailable in this audit.

```sh
python -B -m benchmark_tools.audit_fas_sample_attrition \
  --output NEW_ATTRITION.json --table NEW_ATTRITION.md
python -B -m pytest -q tests/unit/test_audit_fas_sample_attrition.py \
  tests/unit/test_audit_qfo_fas_samples.py \
  tests/unit/test_audit_corrected_fas_samples.py
```

No inference was rerun, no native benchmark score or reference changed, and
no confidence interval is admitted.
