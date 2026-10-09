# Controlled Fragment Representative Trace Protocol

This finite retrospective diagnostic advances publication requirement 4.3.
It uses the completed fragment20_center60_v1 condition and ten unchanged
baseline controls. Aggregate scores have already been exposed; this is not
prospective accuracy validation, a new endpoint or default optimization.
File schemas were inspected before this protocol, but individual error
selections and stage explanations have not been inspected.

## Frozen Inputs

- controlled_fragment_results_20261009_v1/report.json:
  b1480e8eb2ace2aed743f5fc4e515e27ac4ba68b48e644b33989eda3cf8a7ab7
- controlled_fragment_execution_20261009_v1.json:
  af665df1fff1ed454c9cce0286d9bcac831e935ad50b65aba5ee14e76d58edcf
- benchmarks/work/controlled_fragment_observations_20261009_v1/manifest.json:
  8e58830911f60de5aec9ad616e5c02c12229702780e26216aa1371017637666d

Use all seeds 20261101--20261110 and methods in this order: high-sensitivity
OrthoHMM, satellite_v2 phylogenetic OrthoHMM, full OrthoFinder 3.1.5, and its
sequence-only MCL checkpoint. The checkpoint remains diagnostic, not a fourth
fresh inference or an independently timed workflow. Reuse admitted native
prediction adapters and original execution-inventoried artifacts. Bind truth,
flags, input ownership, prediction files and selected stage files to their
retained hashes. Validate the affected scope, not every unrelated old panel.

## Fixed Selection

Canonicalize cross-species gene pairs without changing IDs. Independently
enumerate native baseline and fragment prediction sets with the already
audited converters; their TP/FP/FN counts must reproduce the retained report.
No new F1, confidence interval, inference or scientific admission is created.
Use the original zero/one/two flagged-endpoint strata, including inter-origin
false positives. Apply these seven categories separately within each method
and stratum across all ten seeds:

- fragment_fn: truth minus fragment predictions.
- fragment_fp: fragment predictions minus truth.
- new_fn: truth in baseline predictions but absent from fragment predictions.
- new_fp: false pairs present in fragment but absent from baseline predictions.
- recovered_fn: truth present in fragment but absent from baseline predictions.
- removed_fp: false pairs present in baseline but absent from fragment predictions.
- retained_tp: truth predicted in both arms, a non-error comparison control.

For each of the 4 x 3 x 7 = 84 planned bins, report its complete eligible
pair-record count (seeds are distinct records). Select exactly one record
when nonempty: smallest SHA256 of UTF-8
`method:category:seed:gene_a:gene_b`, breaking ties by (seed, gene_a, gene_b).
Empty bins remain unavailable with count zero and no substituted example.
Selections may overlap across categories; report repeated identities and do
not treat these representatives as independent findings or prevalence estimates.
Do not switch rules after inspecting identities, choose only successes, or
choose cases based on a preferred pipeline explanation.

## Retained Stage Trace

For each selected method/pair, compare baseline and fragment observations.
Also retain prediction presence in all four methods in both arms, allowing
the full OrthoFinder comparator to agree or disagree without being ground truth.

For OrthoHMM inspect directed significant-hit records/scores, final direct
graph edges, graph connectivity and seed/group membership. For satellite_v2,
also inspect candidate co-membership, merge/seed sidecars, final RootHOG
co-membership and the observed pair LCA in saved reconciliation nodes. Reuse
the existing validated graph, partition, node and constraint helpers. Identify
single-copy/unambiguous bypasses explicitly; never invent a reconciliation
node for them. Validate the selected native pair against the recorded event
rule and applicable final root-membership constraint, preserving all observed
stage transitions. Trace the existing MCL group/native pair transition for
OrthoFinder; absence of a matched OrthoFinder search/node adapter is unavailable
stage evidence, not evidence that no hit or duplication exists.

No search, alignment, gene/species-tree inference, clustering or reconciliation
is rerun. No existing predictions or histories are modified. Save selections
before inspecting individual stage explanations. Execute tested deterministic
postprocessing at fresh paths only, after focused commit/push milestones.
Read back selections/counts and stage observations independently, then integrate
the bounded findings and remaining unexplained transitions into the manuscript.

## Interpretation Boundaries

Graph connectivity, co-membership and significant search hits are not orthology
predictions. Missing hits do not distinguish prefilter rejection from weak
scoring. Stage location is an observed exclusion/retention, not a causal
biological explanation or proof that an inferred tree is correct. Do not call
an unchanged score or absent merge a general HMM benefit. Initial HMM search
remains on; there is no HMM-off causal control here. This condition does not
supply natural fragment truth, validated full-domain architecture, independent
biology or publication readiness. The changed full OrthoHMM inventory and
shared-host timing limitations remain; this read-only diagnostic supplies no
new isolated timing comparison or paired baseline resource measurement.
