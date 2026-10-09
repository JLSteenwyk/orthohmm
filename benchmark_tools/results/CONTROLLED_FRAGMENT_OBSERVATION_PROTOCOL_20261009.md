# Prospective Controlled Fragment Observation Test

## Purpose And Exposure

Publication goal 4.3 asks for fragment-related error analysis. Existing missing20
conditions remove genes, and existing variable-length simulations change family
lengths; neither establishes a known truncated observation. This bounded test
addresses that gap by truncating simulated full observations with known parent
sequences. It is not a new method, independent biological validation, a claim
of natural fragment truth, or an indel/domain-architecture simulation.

Baseline method outcomes and their interpretation are already inspected and
development-exposed. No controlled-fragment sequence has yet been generated
or scored. Freeze this protocol and tested adapters before new fragment
inference. Do not change method parameters, seeds, severity, sampled genes,
endpoints or comparisons after seeing fragment outcomes.

## Retained Baselines

Use ALL ten baseline datasets from the existing variable_length_v2 panel,
seeds 20261101 through 20261110 inclusive. No selection by outcome or success
is introduced: the retained result has all four methods admitted for each
baseline seed. Reuse those controls rather than rerunning a completed panel.
Current discovery establishes file presence and report identities; it is not
yet a fresh validation of every raw input, native prediction or runtime.

| Artifact | Relative Result Path | SHA256 |
| --- | --- | --- |
| Retained results | simulation_variable_native_results_20260916.json | cc99fc31c3433809098212d3c8dc12f4ad829b4cfeb66dabcb13aede4e95258f |
| Corrected native method manifest | publication_variable_native_methods_20260916.json | bf677728cf9312edd9755d9ac1ceefbce3abb8c3848d43413631dc00ca00eb2f |
| Generation manifest | publication_variable_simulation_manifest_20260916.json | 806aa1e5f6976c323ff2f6641dd88e25264e7417749e7d77565294666aee776b |
| Reused original comparator manifest | publication_variable_methods_20260916.json | f68b0caf508fde8f40d311b5d303f7569a118bd15e8cf960d227a7d1cfee5984 |
| Native runtime | publication_native_runtime_20260916.json | aebea83807356b02307473506fa30c2dbd2c511d7ba75a0655eb12180a474d74 |

Before consumption, validate each baseline input/truth/prediction against its
actual retained execution/admission bindings and require identical gene IDs,
species owners and truth across methods. Reject missing, changed or unsupported
records; do not silently substitute an output or infer a missing score.

## One Fixed Observation Condition

Name the condition `fragment20_center60_v1`. Within each seed rank every
extant gene ID by SHA-256 of the UTF-8 string
`seed:fragment20_center60_v1:gene_id`, breaking hash ties by gene ID. Select
the first floor(N/5) genes. Selection does not inspect sequence content,
family/truth labels, baseline predictions or accuracy. All ten seeds remain.

For a selected sequence of length L, retain floor(3*L/5) residues from its
center: start=floor((L-retained_length)/2), stop=start+retained_length,
using zero-based half-open coordinates. Preserve every unselected sequence
exactly. Require 0 < retained_length < L; unsupported parent lengths are
reported, not rescued with a post-outcome minimum. Keep every gene/species
ID, species, family membership and event-derived ortholog truth unchanged.
Record original/observed sequence hashes, coordinates and fragment status
for every gene; independently verify the exact substrings and unchanged genes.

Do not regenerate simulator histories, remove genes, alter genealogy, select
an advantageous domain or relax candidate/score thresholds. These are defined
synthetic observation truncations, not biological losses or domain evolution.

## Frozen Methods And Execution

Use the existing frozen OrthoHMM core revision
`7f3a9e40dd7e79f842cc2c11fb8b548f9a802806`, corrected native/profile runtime,
high_sensitivity and satellite_v2 configurations, and full OrthoFinder 3.1.5.
Retain all scientific arguments from the admitted baseline manifest. Change
only input/output/metrics paths for the genuinely unrun fragment identities.
Extract the sequence-only MCL checkpoint from the same full OrthoFinder run
as a diagnostic, with no independent timing claim.

There are 30 new inference identities: ten seeds times three native methods.
Within each seed use the existing high-sensitivity, satellite_v2, full
OrthoFinder order; process seeds ascending and run methods sequentially.
Use the existing scheduler with feasible enforced limits covering the frozen
four-worker/four-threads settings. Test and record exact entrypoint/runtime,
resolved dependency, command and accounting requirements before launch.
Prospective environment compatibility recovery may preserve frozen scientific
settings, but must be tested and explicitly labeled rather than asserting
historical runtime equality. Do not change OS/services or unrelated processes.

Shared Threadripper contention is accepted and disclosed, not a gate. New
resource observations are descriptive and separate from cached baseline
costs; no paired runtime, isolated speed ranking or simulator-cost claim.
No duplicate launches or automatic inference retry. Native failures stay
missing and affect only their dependent conversion/admission/comparisons.

## Scoring And Strata

Reuse the validated cross-species event-derived truth endpoint and pair
conversion semantics. Score all predicted cross-species pairs in the complete
retained gene universe, including false positives between origin families.
Reject unknown IDs and inconsistent species ownership; explicitly canonicalize
and deduplicate orientation. Full orthogroup and native resolved-pair
predictions remain distinguished; the MCL checkpoint is diagnostic only.

Primary summaries are seed-level pair micro F1, precision and recall for all
four methods on the fragment observations. Report prediction volume, gene/pair
coverage and every failure. Verify counts and ratios independently.

Additionally classify every truth and prediction pair by whether exactly zero,
one or two endpoints have the prespecified fragment flag. These three strata
partition all pairs, including unrelated-family false positives. Report all
TP/FP/FN counts, eligible truth and metrics by stratum for baseline and fragment
conditions under the SAME flags. Do not condition on successful recovery or
drop missing satellites. Zero denominators remain explicitly undefined, not
fabricated perfect scores. Stratum differences are descriptive; no subgroup
confidence intervals or outcome-chosen subdivisions.

## Paired Uncertainty

Five planned comparisons, each for F1, precision and recall (15 endpoints):

1. Fragment minus cached baseline for high-sensitivity OrthoHMM.
2. Fragment minus cached baseline for satellite_v2 OrthoHMM.
3. Fragment minus cached baseline for full OrthoFinder.
4. High-sensitivity OrthoHMM minus full OrthoFinder on fragment observations.
5. Satellite_v2 OrthoHMM minus full OrthoFinder on fragment observations.

Recompute the mean paired seed-level metric difference in each of 20,000
whole-seed bootstrap replicates, PCG64 seed20261011, NumPy linear percentile
quantiles. Use nominal 0.025/0.975 and separate fixed-15-endpoint Bonferroni
quantiles 0.05/30 and 1-0.05/30. Reset the RNG per comparison so identical
eligible seed lists share draws across metrics. Use only successful paired
seeds with defined metric values, explicitly retaining excluded IDs/reasons;
do not shrink the planned multiplicity inventory. No available seeds means
unavailable estimates, not zero. Do not subtract means over different seed sets.

These intervals are conditional approximations over the retained independent
simulation master seeds. Ten seeds and finite tails limit coverage, especially
with failures. They are not exact simultaneous, pair-IID, selection-adjusted,
independent biological or generalization guarantees. Existing 14-endpoint
simulation and QfO interval families remain unchanged and separate.

## Deliverables And Boundaries

Commit/push the protocol and tested transformation, new identity/launch adapter,
scorer and summary before production use at fresh output paths. Generate a
single complete manifest for all ten transformed datasets and 30 identities.
Retain baseline pins, parent input/truth, chosen IDs/coordinates and actual
new runtime/command/terminal outcomes. Independently check transformed truth,
native admission, counts, seed comparisons and all tables. Integrate negative,
neutral and failed results and their limits in the manuscript/claim register.

This completes one prespecified controlled-fragment question, not a new tuning
campaign. It does not resolve natural fragment truth, complete domain histories,
all QfO uncertainty or general OrthoHMM superiority. Do not rebuild archives,
repeat the completed simulation panels or alter previously bound evidence.
The full publication goal remains active beyond this diagnostic.
