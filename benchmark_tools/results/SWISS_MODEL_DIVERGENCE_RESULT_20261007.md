# SwissTrees Model-Based Divergence: Bounded Result

## Design And Verification

This advances publication goal4.3 with an explicit model-based descriptor,
rather than calling composition or pair identity a calibrated distance.
[Prospective protocol](SWISS_MODEL_DIVERGENCE_PROTOCOL_20261007.md) was pushed
in24c6e8ae before inference; tested feature/batch/edge-reader sources were
pushed8d075769 before the only selected job23932. All18 retained MAFFT reference
alignments/563 canonical members were reused without extraction, realignment,
masking or outcome-based filtering. Fixed IQ-TREE3.0.1 WAG+G4, seed20261007,
one thread and identical-tip retention were applied to every family.

Job23932 completed0:0 in00:50:38 with2allocated CPUs/8GiB. All18families succeeded
without timeout/retry. Slurm supplied no usable TotalCPU or MaxRSS for this
auxiliary job: these are unavailable, not zero resource use. Its elapsed time
is descriptive feature construction on a shared host, not OrthoHMM inference
cost or a matched, isolated speed measurement. Native23902/review23910 were
left unchanged; cell7's failed timing remains ineligible.

Feature report ceb4ca36/94010bytes is retained inside the small
[input/tree/log evidence archive](swiss_model_divergence_evidence_23932_v1.tar.gz).
The [independent edge-split readback](swiss_model_divergence_readback_23932_v1.json)
45cb12b3/69384bytes verified all18families, all563tips and ALL10765 unordered
tip-pair distances by summing separating edges, not the exporter's path-distance
implementation.190current direct bindings were checked, not the historical
194-file admission or original proteomes. Both implementations share Bio.Phylo
parsing; this is independent arithmetic, not independent tree/model inference.

Family median pair distances range0.463992 to4.129497 expected amino-acid
substitutions/site under the fixed model. The prospectively defined ordinary
median across family medians is approximately2.24582463005. The <=median tie
rule gives9lower/9higher families. No alternate cutoff or bootstrap was tried.
Floating summation differences at about1e-15 are below the reader's1e-10
absolute tolerance and did not change any bin. These are estimated branch-path
distances, not dated evolutionary time or known ancestral history.

Separate projection/rational-reader sources88208814 were committed/pushed before
any selected subgroup score. The [complete table](swiss_model_divergence_strata_20261007_v1/TABLE.md)
and machine-readable reportddfa2775 retain9score rows,6conditional differences,
all54original family-count rows including TN and the same native identities.
[Independent exact-arithmetic readback](swiss_model_divergence_strata_readback_20261007_v1.json)
186c0151 verified every score/difference and human/TSV row. The statistic is
macro-family precision/recall followed by harmonic F1 with the original
TP/2+1, FP/2+1, FN/2+1 prior; neither pooled-pair F1 nor mean-family F1.

## Descriptive Effects

All three cells use initial HMM search and disable downstream profile
refinement. R_at_P0_C0 is the reconciliation effect conditional on candidate
expansion off. C_at_P0_R0 is the expansion effect conditional on reconciliation
off. These three cells do not identify their interaction.

| Distance Bin | Families | Contrast | Delta F1 (pp) | Delta Precision (pp) | Delta Recall (pp) |
| --- | ---: | --- | ---: | ---: | ---: |
| All | 18 | Reconciliation | +10.039 | +30.521 | -6.533 |
| Lower/equal median | 9 | Reconciliation | +11.945 | +32.417 | -7.950 |
| Higher median | 9 | Reconciliation | +8.427 | +28.625 | -5.117 |
| All | 18 | Candidate expansion | -0.347 | -3.523 | +4.375 |
| Lower/equal median | 9 | Candidate expansion | -1.245 | -4.547 | +5.233 |
| Higher median | 9 | Candidate expansion | +0.324 | -2.500 | +3.517 |

Expansion increases recall and decreases precision in both bins; the small
higher-distance F1 increase coexists with a lower-distance decrease and an
overall decrease. Reconciliation increases precision/F1 while reducing recall
in both bins. No subgroup intervals, significance or comparative causal claim
is supported. These patterns do not authorize a distance-dependent default.

## Retention And Limits

Archive639031bytes/SHA256
eeed6e6786f5720ee8c19fb6123fad64c38fdb25c33afe57f5e61bcae4a6113d
contains the18selected alignments and complete new feature attempt outputs,
logs, pair TSV and receipts. It is an evidence supplement, not a replacement
for rc5 or a hermetic release of every project dependency. Existing manuscript,
archive, scientific source and benchmark bytes remain unchanged.

One feature readback took0.47s/57456KiB RSS; one score projection0.05s/13824KiB;
one rational readback0.06s/15360KiB. All exited0 with0swaps. These are auxiliary
postprocessing observations, not inferred tool speed or repaired timing.
Fixtures used wholly invented counts/features, plus a small installed-binary
identical-tip smoke.29feature and18projection tests passed before their selected
executions;32goal/shared-policy tests also passed. No selected attempt failed.

Families are development-exposed, not independent confirmation. Fixed-model
fit, alignment quality and heuristic topology uncertainty are not validated by
execution success. Distances include speciation/paralogy, within-species pairs,
sampling, domains/length, composition and gap/alignment effects. No calibrated
clock, ancestral-history truth, causal explanation, new confidence interval,
default promotion, tuning, accuracy/resource admission or general superiority
is claimed. The full publication goal, other QfO uncertainty limitations and
original TreeFam-artifact gap remain open; this bounded supplement is progress,
not publication readiness.
