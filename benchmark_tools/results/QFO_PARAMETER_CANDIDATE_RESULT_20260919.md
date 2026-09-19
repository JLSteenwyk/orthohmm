# Corrected QfO Candidate Neighborhood Validation

Preparation job 21927 completed with exit 0:0 in 7:40 on bizon with 2 CPUs
and 64 GiB. Independent admission job 21929 completed with exit 0:0 in 1:46,
using the detached executor at `6a0258c15595775752bff5528d7ba678a66725df`.
All five candidate arms passed provenance, parameter, content and merge
reconstruction checks. No orthology accuracy was evaluated.

| Arm | Seed groups | Candidate groups | Reconstructed merges | Input proteins |
| --- | ---: | ---: | ---: | ---: |
| control | 391908 | 351739 | 40169 | 984137 |
| norm_low | 391908 | 351199 | 40709 | 984137 |
| norm_high | 391908 | 352300 | 39608 | 984137 |
| margin_low | 391908 | 345208 | 46700 | 984137 |
| margin_high | 391908 | 355871 | 36037 | 984137 |

These are intermediate candidate-family counts, not final orthogroup counts,
ortholog-pair counts or reference-relative precision/recall. More merges do
not by themselves establish greater sensitivity or better inference.

The unchanged control's candidate partition and complete merge trace match
the admitted corrected baseline byte-for-byte. All arms retain the complete
input universe, use the same profile-refined seeds, and change only the
prespecified scalar. Applied parameters are recorded separately from the
unchanged wrapper profile report. The numeric checkpoint was independently
rechecked: 90,687,327 directed hits, 984,137 genes, 78 species, and no
nonpositive or nonfinite scores. Retained runtime/source/input hashes were
checked before and after preparation; the admission verified those records
and independently reconstructed each candidate partition from its merge trace.
This does not independently recompute search evidence for each merge.

## Retained Evidence

- `qfo_parameter_candidates_prepared_21927.json`, SHA-256
  `e1511d1d8ae6d88f1914ffcba0d09e22ced8253650bee6cedcaa38b8c26f21d2`.
- `qfo_parameter_candidate_admission_21929.json`, SHA-256
  `77759ee0f733242d6d4c368f44c2c1fd84ef3b1ebbb679494acce7fab68f0141`.

Both retained JSON files are byte copies of the completed work artifacts.
Their absolute paths preserve original provenance, not portable acquisition
locations. Large partitions, traces and hit arrays remain local.

GNU time recorded 459.60 seconds wall time, 425.88 seconds user CPU,
33.60 seconds system CPU and 11,192,840 KiB maximum RSS for preparation.
These shared-host cached-preparation measurements exclude initial search,
phylogeny and scoring, and are not controlled end-to-end resource evidence.
The raw record is `benchmarks/work/qfo_parameter_candidates_21927.time`.
The pre-Python failed submission 21926 remains documented in
`QFO_PARAMETER_CANDIDATE_SUBMISSION_21927.md`; it was not erased by success.

## Next Gates

Run each changed candidate partition through its own inferred phylogeny,
with only exact-membership/input/tool-validated raw-tree checkpoint reuse.
Independently admit native groups and cross-species pairs before official
QfO scoring. The full baseline output remains the comparator for all variants.
The two CPM variants still require complete grouping/profile replay from
the fixed hits; existing checked clustering workers bind resolution 0.1,
so explicit tested parameter support is required before those variants run.

The fixed six-variant panel and 18-endpoint SwissTrees multiplicity remain
unchanged. No default was selected, no method superiority was inferred,
and neither the robustness panel nor the publication package is complete.
