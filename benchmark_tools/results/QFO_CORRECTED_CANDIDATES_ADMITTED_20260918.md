# Corrected QfO Candidate Admission

Preparation job 21758 completed successfully in 3:39; independent admission
21759 completed successfully in 44 seconds. All four candidate arms retain
the complete 984,137-protein universe. Admission checks 444 provenance
records and independently verifies recorded merge/partition consistency.

| Profile expansion | Candidate expansion | Seed groups | Candidate groups | Recorded merges |
| --- | --- | ---: | ---: | ---: |
| Off | Off | 394,328 | 394,328 | 0 |
| Off | On | 394,328 | 353,638 | 40,690 |
| On | Off | 391,908 | 391,908 | 0 |
| On | On | 391,908 | 351,739 | 40,169 |

Profile-off retains the initial HMM search; it is not an HMM-free control.
These counts are preparation outcomes, not accuracy endpoints or evidence
that merges are biologically correct. The audit does not independently
recompute the underlying search-support values.

The [unmodified admission receipt](qfo_corrected_candidate_admission_21759.json)
is copied from `benchmarks/work/qfo_corrected_candidate_admission_20260918.json`.
Both 128,164-byte files have SHA-256
`7ae64b1ccd398a7c61c6011f89582c37be91fb9b3954edf85f0cf7d956890f2d`.
The prepared manifest SHA-256 is
`d8385c50426e690afd6d32f3c5302e678de6977c9841451d013201b0f75b564a`.

The first reconciliation job (21760_0) and four non-reconciled pair conversion
jobs (21765, 21767, 21769, 21771) started after admission. Remaining
reconciliations are serialized by the frozen array concurrency limit.
Native reconciliation, conversion, scoring and independent assessment still
require their own admission. No score or uncertainty result follows from
this preparation milestone. Shared-host preparation elapsed time is not a
dedicated end-to-end resource comparison.
