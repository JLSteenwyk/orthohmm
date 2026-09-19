# Corrected Comparator Uncertainty Execution

This implementation note is frozen before calculating corrected whole-panel
intervals. Corrected point estimates and composition-stratified intervals
are already known; this is not a new independent or unexposed test.

The corrected-release protocol explicitly retains the original SwissTrees
comparator design: 18 families, 100,000 shared PCG64 draws, seed 20260920,
eight contrasts, F1/precision/recall and 24-endpoint Bonferroni correction.
Neither protocol is changed. Protocol hashes are checked by the driver:

- `QFO_CORRECTED_RELEASE_PROTOCOL_20260917.md`:
  `b3603bc9b51ce1708f49a0ee816eb02d3cb52c66fb00f07eeb0d66f840a03e1e`.
- `QFO_SWISS_COMPARATOR_UNCERTAINTY_PROTOCOL_20260917.md`:
  `a1ba79268b732824e0df9c56dc1b14dce7b9f33c955523e8b87133c5440161b8`.

The initial comparison inventory is
`qfo_corrected_comparison_20260919_v5/manifest.json`, SHA-256
`7256920214d2eeda9cf596c7c39212992c7f3fcc36dd61e1848fa44d5725e01e`.
It supplies six independent score admissions and explicitly unadmitted
FastOMA/OrthoMCL rows. These missing methods are not imputed, substituted
from historical inputs or silently dropped. The corresponding contrasts
have null intervals, and the multiplicity denominator remains 24.
The complete panel remains incomplete until every planned contrast is
supported by admitted corrected evidence. Missing is not a failure claim.

The driver reconstructs corrected raw counts through the existing admission
bindings, verifies the exact publication rows and native score arithmetic,
requires the same truth identities and family memberships, then computes
macro precision/recall followed by harmonic F1 within each resample.
Historical raw data supply reference labels only. Native phylogenetic pairs,
high-sensitivity group cliques, the OrthoFinder MCL diagnostic and FastOMA's
supplied-tree design remain explicitly distinguished.

All results, including unfavorable and null contrasts, must be retained.
No parameter changes, outcome-selected runs, release-effect tests or
general-superiority claims are authorized. Uncertainty applies only to the
estimated corrected SwissTrees contrasts, not other endpoints or the
secondary mean. Independent numerical reproduction and table/figure
generation remain required after the scheduled run.
