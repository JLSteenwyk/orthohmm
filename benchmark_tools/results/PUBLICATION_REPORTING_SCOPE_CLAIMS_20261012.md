# Publication Reporting Scope And Claim To Evidence

This supplements the [original checklist](PUBLICATION_CLAIMS_20260916.md),
[four-cell checklist](PUBLICATION_FOUR_CELL_CLAIMS_20261007.md) and
[fragment/statistical checklist](PUBLICATION_CURRENT_EVIDENCE_ADDENDUM_20261010.md).
Their dated bytes remain unchanged. Statements that FAS uncertainty is wholly
unavailable predate the conditional sampling analyses below. This October 12
successor integrates the user-approved VGNC/GO/EC, TreeFam and FAS/secondary
reporting amendments while preserving the
[October 10 FAS checklist](PUBLICATION_NATIVE_FAS_CLAIMS_20261010.md).
Biological generalization uncertainty and unconditional historical interval
admission remain unresolved, not completion gates for the explicitly bounded
reporting endpoints. No missing score or failed timing is repaired.

The target is expected repeated native post-attrition FAS under a fixed design,
not the observed score, an unseen-family effect, orthology F1 or the six-metric
secondary mean. September comparators and newer native ablations are separate
snapshots; their ranges must not be combined into a new comparator experiment.

## Claims

| Claim | Evidence | Status And Limit |
| --- | --- | --- |
| Conditional FAS sampling ranges exist for all eight retained September methods | [Frozen reporting plan](NATIVE_FAS_CONDITIONAL_REPORTING_PROTOCOL_20261010.md), [eight ranges and all 28 differences](native_fas_conditional_report_20261010_v1/report.json), [interpretation](NATIVE_FAS_CONDITIONAL_RESULT_20261010.md) | Supported under fixed eligible populations, retained full precomputed means, uniform without-replacement selections with independent strata within a method, and fixed numeric-return status/value. A single joint rectangle supplies all 28 contrasts without cross-method independence. Observed Z remains descriptive; unknown missing-population return count/mean precludes an exact expected-ratio point |
| Both retained OrthoHMM modes exceed full OrthoFinder on this conditional FAS target | [All retained contrasts](native_fas_conditional_report_20261010_v1/summary.md), [separate execution and 80-digit readback](native_fas_conditional_execution_20261010_v1.json) | Supported conditionally: HS minus full OF range [0.076417827, 0.091856675]; PHY minus full OF [0.067014510, 0.075477184]. Adverse contrasts are retained: Proteinortho exceeds both OrthoHMM modes, and HS exceeds PHY. SonicParanoid minus OrthoMCL includes zero. This is not overall accuracy, superiority across benchmarks or causal phylogeny evidence |
| The newer four native cells have conditional FAS ranges despite unknown means in both strata | [Frozen two-sample derivation](NATIVE_FACTORIAL_FAS_SAMPLING_PROTOCOL_20261010.md), [four ranges and all six contrasts](native_factorial_fas_sampling_20261010_v1/report.json), [result and limits](NATIVE_FACTORIAL_FAS_SAMPLING_RESULT_20261010.md) | Supported for P0/C0/R0, P0/C0/R1, P0/C1/R0 and P1/C0/R1 under the explicit design. Three components per cell receive error 0.05/12; eight-corner projection preserves the random denominator. Initial HMM search stays on. R-off group-clique and R-on native pairs have different semantics. No full precomputed mean from the older snapshot is substituted |
| Reconciliation or additional profile refinement has a resolved expected-FAS benefit in the newer cells | [All six native-cell differences](native_factorial_fas_sampling_20261010_v1/summary.md), [independent recurrence and actual execution](native_factorial_fas_sampling_execution_20261010_v1.json) | Not established. P0/C0/R1 minus P0/C0/R0 is [-0.006273752, 0.026395633], despite observed increase 0.010062087. P1/C0/R1 minus P0/C0/R1 is [-0.008154440, 0.008169344]. Four of six ranges include zero, not equivalence. The two excluding-zero ranges compare R-on with candidate-expanded R-off, changing more than one component; they do not isolate a causal reconciliation effect |
| Numerical and finite-design checks prove unconditional historical or biological coverage | [Finite native context result](NATIVE_FAS_CONTEXT_RESULT_20261010.md), [finite design validation](NATIVE_FAS_SAMPLING_VALIDATION_RESULT_20261010.md), [newer native result](NATIVE_FACTORIAL_FAS_SAMPLING_RESULT_20261010.md) | Not established. Context checks cover their selected pairs only; historical RNG state and complete pair context are unavailable. September historical parser/database identities are unestablished; four newer logs are newly pinned, not historically admission-bound. Numeric readback and finite valid-model/control tests do not prove these assumptions for every historical execution. No pair-IID or missing-at-random argument is used |
| The actual figure distinguishes observed scores from expected-ratio ranges | [Caption and scope](NATIVE_FACTORIAL_FAS_FIGURE_20261010.md), [plotted values and assets](native_factorial_fas_sampling_figure_20261010_v1/figure.json), [actual rendering/readback and visual record](native_factorial_fas_figure_execution_20261010_v1.json) | Supported for all four cells and six contrasts. Open circles show offset descriptive native Z/differences, not fitted expected points or centers of biological error bars. Zero-overlap, missing cells and failed timing remain explicit |
| The October 10 manuscript review includes both FAS analyses and the actual figure | [Current manuscript](publication_native_fas_figure_manuscript_20261010_v1.md), [34-page PDF](publication_native_fas_figure_manuscript_20261010_v1_print/document.pdf), [actual execution and selected-page inspection](publication_native_fas_review_execution_20261010_v1.json) | Supported: 41 focused tests pass; 19 citation IDs and 153 unique tracked local assets resolve. All-page text/image bounds scan finds zero violations. Ten pages 11-17/31-33 were actually inspected, including both FAS tables, the full 23-row profile table and FAS image/caption. This is scoped review, not all-page visual certification, independent science or whole-study native restoration |
| These outcomes establish publication readiness or resolve other QfO uncertainty | [Current manuscript limitations](publication_native_fas_figure_manuscript_20261010_v1.md), [original full requirement register](PUBLICATION_REQUIREMENT_AUDIT_20261004.md), [retained VGNC methodological decision](VGNC_UNCERTAINTY_METHOD_REVIEW_20260927.md) | Not established. Biological/family generalization uncertainty, VGNC/TreeFam/GO/EC and secondary-mean uncertainty remain unresolved. Original TreeFam family resources remain unavailable; three native score cells and failed/ineligible timing stay missing. Negative F1 comparisons remain unchanged. No general OrthoFinder superiority, isolated efficiency ranking, new default or archival deposition follows |

## Approved Reporting Claims

| Claim | Evidence | Status And Limit |
| --- | --- | --- |
| Native VGNC precision/recall and GO/EC bars can accompany descriptive benchmark results | [Revised manuscript and native-bar tables](publication_reporting_scope_manuscript_20261012_v1.md), [original corrected register](qfo_corrected_comparison_20260926_v7/manifest.json), [retained GO/EC runtime audit](QFO_GO_EC_ARITHMETIC_AUDIT_20260917.md) | Supported as conventional reporting: all 24 aggregate hashes and metric values checked against retained admissions. VGNC uses 1.96*sqrt(p*(1-p)/n); GO/EC's native StdErr is a Student-t 95% mean half-width, not one SEM. Shared-protein/family dependence is unaccounted for; no paired significance or propagated F1 interval |
| TreeFam scores can be retained as supporting evidence without source-tree recovery | [TreeFam search scope](TREEFAM_REMAINING_CONTEXT_RESULT_20260930.md), [corrected scores](qfo_corrected_comparison_20260926_v7/scores.md), [revised limitations](publication_reporting_scope_manuscript_20261012_v1.md) | Descriptive supporting evidence under explicit user authorization. Missing original trees/mapping and family-level intervals are documented limitations, not completion gates. No reconstructed family labels, TreeFam significance or global source-unavailability claim |
| Existing FAS ranges can be reported without biological/generalization intervals | [Eight-method conditional ranges](native_fas_conditional_report_20261010_v1/report.json), [newer four-cell ranges](native_factorial_fas_sampling_20261010_v1/report.json), [revised caption and Methods](publication_reporting_scope_manuscript_20261012_v1.md) | Supported only for each stated conditional expected-statistic target, assumptions and snapshot. Observed Z is separate. Historical bindings and sampling/return assumptions remain limitations; not intervals centered on observed scores or general biological superiority |
| The six-metric mean can be shown as a descriptive secondary summary | [Unchanged aggregation](qfo_corrected_comparison_20260926_v7/manifest.json), [revised Methods](publication_reporting_scope_manuscript_20261012_v1.md) | Descriptive only, alongside all endpoints and coverage. Not an official QfO accuracy statistic, significance test, headline superiority or basis for choosing a winner. No aggregate interval or metric-independence assumption |
| The revised manuscript has an actual render and scoped visual review | [37-page PDF](publication_reporting_scope_manuscript_20261012_v1_print/document.pdf), [actual review](PUBLICATION_REPORTING_SCOPE_REVIEW_20261012.md) | Supported: 56 focused tests passed, 19 citation IDs and 156 tracked direct targets resolve, all-page bounds scan finds zero violations, and 13 selected pages were inspected. Not all-page visual certification, new science or a whole-study restoration |
| The amendments alone prove full publication completion | [Amended goal](../PUBLICATION_GOAL_CURRENT.txt), [current 43-bullet reconciliation](PUBLICATION_REQUIREMENT_RECONCILIATION_20261012.md) | Not established by scope edits. Nonblocking gaps do not become solved science. Independent-validation, ablation, robustness, resource, biological-application and reproducible-package evidence are assessed at their actual scope. Historical readiness/uncertainty flags do not override fresh approved reporting scope |

## Scope Of Delivery

This is a post-results, explicitly authorized reporting-scope revision, not
a newly prespecified experiment or a new estimator. Numerical FAS ranges,
native scores, failed outcomes, paired analyses and scientific defaults remain
unchanged. The October 10 review evidence above applies to its dated manuscript,
not a certification of the October 12 successor. New review artifacts are
linked from the reproduction guide when actually completed.

The native-bar readback checks admission-pinned files, axes, participant
consistency and metric arithmetic, not complete transitive identity or
biological interval validity. An initial null-participant assertion stopped
before output; the manifest leaves that field null for five comparators.
Its corrected handling retains each aggregate's actual participant ID and
admission-bound source, checks all three agree, and matches all values/counts
to the existing manifest. No identity normalization or endpoint changed.
