# Four Cell Native QfO Claim To Evidence Addendum

Four admitted native QfO configurations now have checked scores, conditional
SwissTrees uncertainty, a complete profile-pair explanation and reviewed main
text. The evidence supports endpoint-specific trade-offs, not a general benefit
of downstream profile refinement or superiority over full OrthoFinder.
This addendum supplements the [earlier error-stratum checklist](PUBLICATION_CLAIMS_20261007.md)
and [full original claim checklist](PUBLICATION_CLAIMS_20260916.md); it does not
replace them or audit completion of every publication requirement.

## Scientific Claims

| Claim | Evidence | Status And Limit |
| --- | --- | --- |
| Four native QfO configurations have admitted accuracy | [Exact score snapshot](native_qfo_scientific_scores_20261007_v3/report.json), [generation and content checks](FOUR_CELL_MAIN_TEXT_RESULT_20261007.md) | Supported for P0/C0/R0, P0/C0/R1, P0/C1/R0 and P1/C0/R1 across all six endpoints. Three score cells remain unavailable. P is downstream profile refinement, C candidate expansion and R reconciliation; initial sensitive HMM search is on throughout. Native and group-clique output semantics differ; this is not an initial-HMM-off control |
| Prediction coverage is comparable to an accuracy or reference-recall estimate | [Coverage and output semantics](native_qfo_scientific_scores_20261007_v3/report.json), [figure and content readback](native_qfo_four_cell_figure_readback_20261007_v1.json) | Not supported: the common denominator is all 984,137 input accessions. Coverage is participation in submitted relations, not scored-reference recall or correctness. All four conversions have zero mapping loss; this does not validate their biological predictions |
| Three native conditional SwissTrees contrasts can reuse the retained uncertainty analysis | [Guarded binding](native_qfo_four_cell_swiss_uncertainty_20261007_v1.json), [independent rational readback](native_qfo_profile_swiss_readback_20261007_v1.json) | Supported for C at P0/R0, R at P0/C0 and P at C0/R1 after exact native family-record matching. All 100,000 original paired-family draws, seed 20260922 and 42-endpoint adjustment reused. Eleven contrasts unavailable. These are conditional effects in 18 development-exposed families, not independent validation, a complete interaction analysis or pair-independent confidence intervals |
| Downstream profile refinement improves SwissTrees F1 with C0/R1 fixed | [Profile contrast and all family differences](native_qfo_four_cell_swiss_uncertainty_20261007_v1.json), [separate count arithmetic](native_qfo_profile_swiss_readback_20261007_v1.json) | Improvement is not established. Point difference is -0.3202 percentage points; adjusted interval [-2.1129, 0.5132] includes zero. CASP improves, GH14 worsens and 16 families tie. No equivalence, default promotion or unseen-family effect follows. Count-based arithmetic and original serialized endpoints remain distinct |
| Three fewer false positives guarantee a higher macro F1 | [Complete changed-pair localization](native_qfo_profile_localization_20261007_v1.json), [family-level contrast](native_qfo_four_cell_swiss_uncertainty_20261007_v1.json) | Contradicted for this observed contrast: three CASP FPs disappear but one GH14 TP is lost, and macro F1 decreases. Equal family weighting and different precision/recall changes prevent interpreting a favorable pair-count balance as a macro-accuracy gain |
| The four changed assessed profile pairs disappear through candidate separation before reconciliation | [All-pair trace](native_qfo_profile_localization_20261007_v1.json), [independent saved-tree and membership check](native_qfo_profile_newick_readback_20261007_v1.json), [mechanistic result](NATIVE_PROFILE_LOCALIZATION_RESULT_20261007.md) | Supported as a software-stage explanation for every changed assessed SwissTrees pair. Six families, 70 leaves and two before-state LCAs checked independently. Before, each pair shares a family with a speciation pair-event; afterward, endpoints occupy different families and no shared-family LCA exists. This is not a new duplication exclusion, true biological history or proof that an individual HMM edge caused separation. Species-tree byte hashes also differ; no isolated tree effect is identified |
| The earlier fixed-bin and model-distance strata estimate the new profile effect | [Earlier complete strata](NATIVE_QFO_THREE_CELL_STRATA_RESULT_20261007.md), [model-distance results](SWISS_MODEL_DIVERGENCE_RESULT_20261007.md), [current manuscript](PUBLICATION_MAIN_TEXT_20261007_v4.md) | Not supported: those strata cover P0 cells. They remain valid at their original scope but cannot be relabeled as the P1 contrast. No new subgroup interval or profile-specific default follows |
| The new results establish superiority over full OrthoFinder or universal HMM benefit | [Retained all-tool comparison](PUBLICATION_MAIN_TEXT_20261007_v4.md), [four-cell evidence and scope](NATIVE_FOUR_CELL_FIGURE_RESULT_20261007.md) | Not established. The new cells contrast OrthoHMM components, not an independently tuned comparator; initial HMM search is retained. Earlier endpoint-specific comparator evidence and development exposure remain unchanged. The six-metric mean is a project-defined secondary summary, not official QfO F1 |

## Reporting Claims

| Claim | Evidence | Status And Limit |
| --- | --- | --- |
| The current four-cell figure retains all supported endpoints and contrasts | [Manifest](native_qfo_four_cell_figure_20261007_v1/manifest.json), [independent readback](native_qfo_four_cell_figure_readback_20261007_v1.json), [actual visual review](NATIVE_FOUR_CELL_FIGURE_VISUAL_REVIEW_20261007.md) | Supported for 24 score values, four coverage rows and nine interval rows. Actual PNG and decoded figure PDF inspected. Separate F1, functional similarity and coverage scopes remain explicit; it is not a complete factorial or timing comparison |
| The v4 manuscript delivers the new evidence without changing older scientific content | [Exact generation receipt](publication_four_cell_main_generation_20261007_v1.json), [actual content checks](FOUR_CELL_MAIN_TEXT_RESULT_20261007.md), [review result](FOUR_CELL_MAIN_REVIEW_RESULT_20261007.md) | Supported for the documented native section replacement and three narrow ancillary clarifications. Remaining scientific body restores exactly to v3. Original method, all-tool tables, P0 strata, citations and historical archives unchanged. No new inference, scoring, bootstrap or admission |
| The actual new manuscript PDF is readable and linked to its checked source | [Asset receipt](publication_main_review_20261007_v4_checked_assets.json), [print receipt](publication_main_review_20261007_v4_checked_print/print.json), [bounds/page records](publication_main_review_20261007_v4_checked_pdf_review/report.json), [all-page visual record](FOUR_CELL_MAIN_VISUAL_REVIEW_20261007.md) | Supported for all 21 actual PDF pages, zero bounds violations and observed absence of clipping/overlap. Nineteen citation IDs and 107 unique tracked local targets checked. First render's bibliography-selection failure remains failed; the separate corrected render uses the actual parent-bound bibliography. Figure linked, not embedded. This does not establish transitive portable reproduction, scientific correctness or journal formatting |
| This update closes the publication goal or repairs failed timing | [Current limitations](PUBLICATION_MAIN_TEXT_20261007_v4.md), [original checklist](PUBLICATION_CLAIMS_20260916.md) | Not established. Missing native score cells, other-endpoint paired uncertainty, original TreeFam resources, independent-family/history/fragment truth limitations and complete executable archival delivery remain explicit. The failed recovered-cell timing stays ineligible; shared-host observations have unknown, potentially tool-dependent contention effects. No isolated efficiency ranking, publication readiness or archival DOI is claimed |

Actual job states belong in the live ledger and scheduler, not this fixed
scientific addendum. Future native outcomes require their original reviews,
scientific admission and separately versioned reporting; a live run never
supplies a missing score. No automatic retry or reinterpretation of a failed
attempt follows from this checklist.

Checklist verification: Pandoc structural parsing finds eight scientific and
four reporting claim rows. All 32 local link occurrences resolve to 22 existing
tracked targets whose current hashes were checked. The new profile F1 effect,
adjusted endpoints and 1/16/1 family outcome counts match the bound contrast
directly. The previous checklist, selected generator and v4 manuscript have
empty scoped diffs. Link/content checks do not independently validate the
science behind every claim or establish full-goal completion.
