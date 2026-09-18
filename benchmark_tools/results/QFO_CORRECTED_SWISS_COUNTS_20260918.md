# Corrected SwissTrees Raw-Count Audit

The auditor pins independently admitted corrected assessments, checks conversion
bindings, and follows the admitted execution record to its inventoried raw
output. Selected evidence is rehashed before and after reading. Full inference
and scoring provenance relies on the prior admission rather than a new rerun.

All 18 families and 10,765 reference relations agree exactly with the pinned
reference anchor in identities, truth labels and represented members. No genes
are shared between reference families. Historical predictions are not reused;
only reference identity and truth labels are compared. This does not establish
independence of families' evolutionary histories.

| Corrected method | Macro precision | Macro recall | Harmonic F1 |
| --- | ---: | ---: | ---: |
| Proteinortho 6.3.6 | 0.952179936 | 0.576414360 | 0.718110999 |
| SonicParanoid 2.0.9 | 0.872278767 | 0.736159586 | 0.798459419 |

Family confusion counts use `raw / 2 + 1`. Precision and recall are averaged
before calculating harmonic F1, not the mean of family F1 values. Reconstructed
family and aggregate native metrics agree within absolute tolerance 5e-8.
Tiny differences from the score table reflect native decimal serialization.

Count reports and SHA-256:

- `qfo_corrected_proteinortho_swiss_counts_20260918.json`:
  `9813f8adeba9f5149b9a82ee29e8090757282f272652c16e60b45f9f4002886b`.
- `qfo_corrected_sonic_swiss_counts_20260918.json`:
  `8c40764fcc2822e5564652937ea843701dfc3311eaa84e5d0ea00c187afa06c1`.

Reproduce with `benchmark_tools/audit_qfo_corrected_swiss.py --admission PATH
--admission-sha256 SHA --baseline
benchmark_tools/results/qfo_swiss_counts_20260917.json --output FRESH_PATH`.
Raw artifacts must exist at their recorded paths. Reports include exact source,
helper and selected evidence identities.

Nine new tests cover fresh counts, changed reference truth/members, inventory,
coverage, native scores and participant identity. The focused suite including
raw-parser and corrected-score export tests passes 35 tests. Both actual
corrected assessments passed the new count audit.

No intervals or release-effect tests are calculated. The corrected protocol
still requires all eight comparators and its prespecified 24-endpoint family;
do not reduce it to the two complete methods. Corrected factorial uncertainty
is separately prespecified across 42 endpoints. This audit provides no
uncertainty for other challenges or the secondary six-metric mean.
