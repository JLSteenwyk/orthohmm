# Simulation Search-Sensitivity Calibration Protocol

## Scope and Rationale

The existing OrthoBench and QfO DIAMOND controls use nominally equal cutoffs,
not established equal sensitivity. Add a separate search-only diagnostic on
the retained validated variable-length Zombi/Pyvolve simulation panel. Preserve
all previous controls. No OrthoBench/QfO/YGOB score or orthology F1 may be used
to choose the new cutoff; no scientific default changes.

Use all seven existing conditions and ten seeds (70 datasets), without dropping
poor-performing conditions. Seeds 20261101-20261105 form the 35-dataset
calibration split; seeds 20261106-20261110 form the 35-dataset reporting split.
The simulations and past method results are development-exposed: this split is
not independent validation. Conditions with the same seed share histories.

## Frozen Endpoints and Selection

Homology truth is membership in the same simulated ancestral family, including
paralogs. Search is directed; exclude self and within-species hits. A family's
eligible denominator is n squared minus the sum of its per-species counts
squared. Retain zero-denominator families in the inventory but exclude them
from recall means. Primary recall is the unweighted family mean per dataset,
then the unweighted dataset mean. Also report pooled counts and per-condition
recall, nonhomolog hit counts and total reported hits. Homology recall is not
orthology recall or F1. Engine processes receive sequences, not family labels.

Fix the installed frozen-source HMM initial-search settings at BLOSUM62,
k=4, maximum 100 candidates per query, band width 64 and E-value 1e-4, with
four workers/one thread each. This is the initial search of the high-sensitivity
configuration, not multi-sequence profile expansion or reconciliation.

Use the previously pinned DIAMOND 2.1.11 binary, separate target-species
databases, very-sensitive mode, BLOSUM62, gap open/extension 11/1,
composition statistics 1, masking 1, unlimited targets and one HSP per pair.
Search once at E-value 1 with four threads and retain raw scores, bit scores
and E-values. Post-filter at this fixed grid:

`1e-100, 1e-80, 1e-60, 1e-40, 1e-30, 1e-20, 1e-15, 1e-10, 1e-8,
1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1`.

Choose one global cutoff minimizing absolute calibration-split primary-recall
difference from HMM. Exact ties select the smaller E-value. Do not use
orthology grouping scores, F1, reporting-split outcomes or computational cost
for cutoff selection. Show the complete grid, not only the selected point.

The reporting-split match gate requires absolute primary-recall difference
at most 0.02 overall and at most 0.05 in every condition. A failed gate remains
failed; do not extend the grid or choose condition-specific cutoffs after
results. Preserve failures as missing, never as zero hits. Complete-panel
calibration and gate claims require every planned dataset for both engines.
No significance claims or confidence intervals are prespecified for this
feasibility diagnostic.

## Execution and Interpretation

Input preparation must verify historical truth hashes, FASTA hashes, complete
input inventories, unique genes, family coverage and eligible denominators
before search. The installed source/runtime and DIAMOND identity are pinned.
Record worker source, exact commands, environment, input/output hashes,
scheduler state and failures. No automatic inference retries. Run separate
from job 22179 without modifying its inputs, executor or allocation. Shared-host
timing is descriptive; post-filtered broad DIAMOND output does not match effort.

Search recall matching on these simulated histories would be bounded evidence,
not real-data sensitivity equivalence. This panel does not itself estimate
downstream orthology effects: a downstream matched-sensitivity experiment needs
a separately frozen protocol after the match gate is evaluated. The existing
real-data controls and their uncertainty remain unchanged regardless of outcome.

## Verified Input Preparation

[Prepared input/endpoint manifest](search_sensitivity_inputs_20260926.json),
SHA-256 `b3b19927821f46f0e79fcc444987a51b88f699e7665bd364fa6639a4239323bc`,
verifies all 70 planned datasets against admitted historical truth identities
and FASTA hashes. The two splits contain 35 datasets each. Input sizes span
500-1,004 genes; eligible ancestral families span 77-100 per dataset. The sum
of directed homology denominators is 382,854 across dataset instances, not
382,854 independent observations. All 13 preparation tests pass.

Preparation has not executed either new search engine, selected a threshold,
evaluated the reporting gate or changed an existing scientific result. Next
implement and pin the search-only worker, preserving its raw per-species hits,
then launch the complete panel under bounded Slurm resources. The protocol and
input validation are prerequisites, not evidence of matched sensitivity.
