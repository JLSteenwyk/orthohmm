# Complete Corrected QfO Factorial

## Validation

Final-cell scoring job 21787 completed with exit 0:0 in 28:48; independent
scoring admission 21788 completed in 12 seconds. Complete-factorial export,
SwissTrees count audit and bootstrap job 21894 completed with exit 0:0.
Fresh executions of the frozen final-cell admission, count auditor and
bootstrap reproduced all three JSON outputs byte-for-byte. No scientific
configuration, reference, endpoint, seed or multiplicity rule changed.

- [Eight-cell score table](qfo_corrected_factorial_complete_20260919/scores/scores.md)
  and [manifest](qfo_corrected_factorial_complete_20260919/scores/manifest.json).
- [All 42 intervals and family wins/ties/losses](qfo_corrected_factorial_complete_20260919/swiss_bootstrap.md).
- [Machine-readable intervals](qfo_corrected_factorial_complete_20260919/swiss_bootstrap.json)
  and [independently audited counts](qfo_corrected_factorial_complete_20260919/swiss_counts.json).
- [Final-cell scoring admission](qfo_corrected_factorial_score_admission_21788.json).
- [All-cell and all-endpoint figure](qfo_corrected_factorial_figures_20260919/qfo_factorial_swiss.png)
  with [rendering provenance](CORRECTED_QFO_FACTORIAL_FIGURE_20260919.md).

Final `p1_c1_r1` point estimates are GO 0.490349, EC 0.965650,
VGNC F1 0.901690, SwissTrees F1 0.833513, TreeFam-A F1 0.614864 and
FAS 0.762993. The project-defined secondary mean is 0.761510, not a
primary endpoint. All 5,959,560 submitted native pairs map without loss.
SwissTrees precision is 0.955177 and recall is 0.739341. Pair counts
measure prediction volume, not protein coverage.

## SwissTrees Effects

The frozen procedure uses 100,000 shared family-bootstrap draws, seed
20260922, 18 reference families and Bonferroni adjustment across 42
endpoints. Each replicate recomputes the harmonic mean of mean family
precision and mean family recall, rather than averaging family F1.
Differences below are raw 0-to-1 units; intervals are adjusted.

| Contrast | F1 difference | Adjusted interval |
| --- | ---: | --- |
| Reconciliation at P0, C1 | +0.150644 | [0.026173, 0.298879] |
| Reconciliation at P1, C1 | +0.147426 | [0.024122, 0.295846] |
| Candidate-expansion by reconciliation interaction at P0 | +0.050254 | [0.003650, 0.124556] |
| Candidate-expansion by reconciliation interaction at P1 | +0.046552 | [0.001825, 0.121851] |

These four of 14 F1 intervals exclude zero; the other ten include zero.
In particular, all four profile-refinement F1 intervals include zero.
All four reconciliation precision intervals are positive and all four
recall intervals negative. Candidate expansion without reconciliation has
positive adjusted recall intervals at both profile settings, but its F1
intervals include zero. Both F1 interactions have 9 positive, 7 tied and
2 negative family effects; these counts are interaction signs, not a
ranking of methods.

This supports an end-to-end interaction on the development-exposed
SwissTrees benchmark, not an isolated mechanistic effect: expansion can
change species-tree estimation, and R changes group-derived clique pairs
to native phylogenetic pairs. P-off retains initial HMM search. Only 18
families are available; exchangeability may be imperfect. The adjustment
does not cover prior development selection, other benchmarks, or the
secondary mean. No TreeFam family intervals are inferred. No defaults are
retuned, and no superiority over corrected full OrthoFinder is established.

## Checksums

| Artifact | SHA-256 |
| --- | --- |
| Final scoring admission | `49b7d837b2ba2928b0676c9974e2ec16db5211f1086c366c7bc11805a2ee751f` |
| Score manifest | `06c6e0e6497af9e1ed6a0e40d0430601f740f1130d6565d1b5f28aab9cb7d14d` |
| SwissTrees counts | `c9d8bf02ef6f287c56d073fa61f165c982caf834a94fc6e44e6589f03e46eba6` |
| SwissTrees bootstrap | `f777dead1294b3810c0479877aade7fb4411c0544696f671226ff198432e9211` |

Original work-directory paths remain in copied manifests to preserve
provenance; historical partial exports remain unchanged. Full corrected
competitor comparison, composition strata, controlled scaling and wider
publication requirements remain unfinished.
