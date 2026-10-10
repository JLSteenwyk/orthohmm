# Current Fragment And Statistical Claim Addendum

The [current manuscript](publication_supplements_manuscript_20261010_v1.md)
supports a bounded synthetic-fragment comparison and observed software-stage
localization, not natural-fragment superiority. Both HMM configurations have
lower fragment F1 than full OrthoFinder in the tested condition. Relocated
postprocessing is reproducible within its recorded scope; the restricted
Poisson validation does not supply native benchmark confidence intervals.

This supplements the [original claim checklist](PUBLICATION_CLAIMS_20260916.md)
and [four-cell native QfO addendum](PUBLICATION_FOUR_CELL_CLAIMS_20261007.md).
Their dated snapshots remain historical. No result, endpoint, default,
uncertainty admission or publication criterion changes here.

## Claim To Evidence

| Claim | Evidence | Status And Limit |
| --- | --- | --- |
| Controlled truncation was tested without prediction-based selection | [Frozen observation protocol](CONTROLLED_FRAGMENT_OBSERVATION_PROTOCOL_20261009.md), [actual results](controlled_fragment_results_20261009_v1/report.json), [interpretation](CONTROLLED_FRAGMENT_RESULT_20261009.md) | Supported for seeded ID-selected center-60% retention in floor-20% of genes across ten simulation master seeds. All 40 fragment method/checkpoint outcomes were admitted. IDs, owners and evolutionary truth were preserved. This is development-exposed synthetic truncation, not validated natural fragment recovery or domain/indel evolution |
| OrthoHMM outperforms full OrthoFinder on these fragments | [All 15 paired comparisons](controlled_fragment_results_20261009_v1/report.json), [complete comparison table](controlled_fragment_results_20261009_v1/comparisons.tsv) | Contradicted within this condition: HS and PHY both have negative F1 differences and negative adjusted interval endpoints versus full OF. The intervals are conditional seed-bootstrap approximations, not exact simultaneous guarantees or independent biological generalization. Both methods retain initial HMM search; no initial-HMM-off control |
| Small fragment-minus-baseline changes establish robustness or equivalence | [All 80 scores and paired comparisons](controlled_fragment_results_20261009_v1/report.json), [scope](CONTROLLED_FRAGMENT_RESULT_20261009.md) | Not established. All three native F1 change intervals include zero. Inclusion of zero is not an equivalence test, absence of fragment sensitivity or permission for outcome-led tuning. The OF checkpoint is a diagnostic extracted from its full run, not an independently timed sequence-only pipeline |
| Representative traces explain observed software decisions | [53 selected cases and empty bins](controlled_fragment_trace_selection_20261009_v1/selection.json), [106 stage observations](controlled_fragment_stage_trace_20261010_v1/report.json), [independent readback](controlled_fragment_stage_readback_20261010_v2.json) | Supported where saved stages are observed: exact hits, graph connectivity, grouping, candidate membership and saved event/root decisions across 44 contexts. Selection is retrospective and minimum-hash within fixed categories; counts do not estimate error prevalence. Absent significant hits do not distinguish prefilter rejection from scoring/threshold rejection. Missing OF search/event adapters remain null |
| The traces prove a causal HMM defect or erroneous biological duplication history | [Complete observations](controlled_fragment_stage_trace_20261010_v1/report.json), [current limitations](publication_supplements_manuscript_20261010_v1.md) | Not established. A saved duplication call describes the pipeline decision, not evolutionary correctness or an individual causal HMM edge. Connectivity is not orthology; before/after endpoint counts alone do not localize a cause |
| Selected fragment raw-stage postprocessing works at relocated paths | [Actual copied readback](controlled_fragment_relocated_readback_20261010_v2.json), [execution including corrected replay](controlled_fragment_relocation_execution_20261010_v2.json), [preserved first failure](controlled_fragment_relocation_execution_20261010_v1.json) | Supported: 1,783 pinned artifacts and 39,428,185 payload bytes reproduce the full 106-observation result using unchanged kernels with no original-path fallback. Python-open auditing is not OS containment. Raw payloads remain local; this is not historical inference reproduction, public data availability or whole-study restoration |
| The exact count-model construction validates native benchmark intervals | [Frozen model](PAIRED_POISSON_F1_PROTOCOL_20261010.md), [12 cells and two failure controls](paired_poisson_f1_validation_20261010_v1.json), [actual execution/readback](paired_poisson_f1_execution_20261010_v1.json) | Not established. All 12 specified perfect-recall Poisson cases pass, but both non-Poisson common-shock controls have zero coverage. The target is the F1 difference at expected count vectors, not expected sample F1 or inference to unseen families. Native VGNC/TreeFam/GO/EC/FAS and secondary-mean uncertainty remain unresolved |
| The current manuscript includes these outcomes without altering prior science | [Generation record](publication_supplements_generation_20261010_v1.json), [actual scoped readback](publication_supplements_execution_20261010_v1.json) | Supported: two insertions, six new links and exact restoration of the complete 101,043-byte parent after removing those insertions. The earlier 30-page PDF renders the earlier text, not this successor. No native interval, score or readiness is newly admitted |
| These updates complete the full publication goal | [Current manuscript limitations](publication_supplements_manuscript_20261010_v1.md), [live progress ledger](PUBLICATION_PROGRESS.md) | Not established. Native uncertainty, unavailable source resources and failed native outcomes remain explicit. No new isolated timing, independent accuracy confirmation, general OrthoFinder superiority, final scientific admission or archival deposition follows |

## Retained F1 Effects

Values below are generated from the full-precision retained JSON, in percentage
points, rather than recomputed from rounded displayed scores. HS is high
sensitivity, PHY is satellite_v2 with phylogeny, and OF is full OrthoFinder.

| F1 Comparison | Difference (pp) | Fixed-15 Adjusted Interval (pp) | Paired Seeds |
| --- | ---: | --- | ---: |
| HS fragment minus baseline | -0.0999 | [-0.3878, 0.0000] | 10 |
| PHY fragment minus baseline | -0.1223 | [-0.4835, 0.0893] | 10 |
| OF fragment minus baseline | -0.0344 | [-0.1237, 0.0000] | 10 |
| HS fragment minus OF fragment | -2.9307 | [-5.6825, -1.0863] | 10 |
| PHY fragment minus OF fragment | -0.5719 | [-1.9815, -0.0087] | 10 |

All 15 F1/precision/recall endpoints were retained, using 20,000 paired
whole-seed draws with fixed-15 endpoint adjustment. This checklist reuses
those outcomes without new draws or scores. The synthetic seed-level units
are not reference-family or native dyadic resampling units, and do not
justify confidence intervals for the unresolved QfO endpoints.
