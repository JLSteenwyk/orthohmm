# Forced Candidate Diagnostic

This follow-up is specified after inspecting OrthoBench development results.
It does not alter publication defaults or constitute independent validation.

Score all 81,466 directed within-reference-family pairs from the existing
70-family trace, including the 48,368 previously excluded by the prefilter.
Use the frozen publication core, CPU scoring, BLOSUM62, band width 64 and
the original full target proteomes for E-value database size. Override only
the candidate arrays supplied to the existing engine. Preserve its profile
construction, normalization and significance threshold (strict E < 1e-4).
Retain the complete query subset per species direction and all 144 directions.
The temporary candidate override must be restored even if scoring raises.

Before scientific interpretation, independently join every directed pair to
the original decision trace. Require identical universes and compare numeric
scores and E-values for previously scored pairs, not just pass/fail status.
Report accepted/rejected counts separately for prior prefilter exclusions,
prior significance rejections and prior acceptances. Preserve disagreements
and investigate them before attributing any rescue to candidate selection.

Reference labels determine which candidates are forced. Therefore rescued
reference pairs are mechanistic evidence only: do not call their fraction
unbiased sensitivity, predict a whole-proteome precision gain, or modify the
production pipeline to inject reference pairs. This does not establish
matched DIAMOND sensitivity, isolate all differences between search engines,
or measure the downstream accuracy of a broader prefilter. Those limitations
remain even if all previously excluded pairs pass significance.

Results are uncomputed at protocol creation. Shared-host runtime is descriptive.
