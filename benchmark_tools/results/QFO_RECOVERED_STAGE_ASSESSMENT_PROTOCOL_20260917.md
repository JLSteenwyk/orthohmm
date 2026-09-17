# Recovered QfO Stage Assessment

Freeze this panel before scoring the independently admitted recovered QfO replay.
This is development-exposed stage analysis, not independent validation or a
replacement for the historical benchmark rows. Admission SHA256:
905c8accae2a0c98e73e028ec30e8a40d5a3d36d3223169afb6bedc01127758e.

## Fixed Inputs and Semantics

Retain allfour output partitions in their frozen order: multipass,
multipass_refined, strict_profiles, strict_profiles_refined. Use the unchanged
78input FASTAs and full976504-gene universe. Recheck admission, partition,
converter and input hashes. Conversion uses the existing audited
`qfo_benchmark/og_to_pairwise.py` cross-species clique expansion, followed by the
existing `qfo_filter_pairs.py` against the2020mapping. Record both raw and retained
pair counts and mapping exclusions. These are cluster-derived pair predictions,
not native reconciled ortholog pairs. No gene or group is selected by a score.

Use unique participants ohmm_checked_v2_0 through ohmm_checked_v2_3 and isolated
output/work namespaces. Never overwrite historical QfO results. Preserve failed
conversion/scoring attempts and never substitute historical scores for these
new partitions. Preparation does not run inference or compute accuracy.

## Assessment and Interpretation

Evaluate allsix established2020challenges: GO, EC, VGNC, SwissTrees, TreeFam-A,
and FAS, using the existing benchmark-webservice and container versions. Pin the
actual pipeline source, reference files, mapping and local image hashes before
assessment. Use short work paths within the established Darwin path-length limit.
Retain individual native metric names, values, denominators and completion state.
Any six-metric mean remains a project-defined secondary summary, not official F1.

Report the fixed four contrasts: strict_profiles versus multipass,
strict_profiles_refined versus multipass_refined, multipass_refined versus
multipass, and strict_profiles_refined versus strict_profiles. Do not choose the
best stage after seeing results. This crosses profile-branch execution and
sequence-based post-clustering refinement; it is not the candidate-expansion by
phylogenetic-reconciliation factorial still required separately.

The profile branch starts from unrefined multipass clusters and adds HMM-derived
edges to the base graph, reclusters, recomputes singleton-assignment edges, and
reclusters again. Therefore its differences include downstream graph/assignment
responses, not only direct HMM score effects. Post-clustering refinement here is
not phylogenetic reconciliation. The generating replay used an explicitly
checked Python-pair graph constructor and isolated clustering-child thread
settings; do not silently identify it with the historical implementation.

No claim of significance follows from scalar metric differences alone. Paired
uncertainty must respect each challenge's reference units and statistic; where
the necessary per-unit evidence is absent, report uncertainty as unavailable
until a validated challenge-specific analysis exists. Dependent predicted pairs
must not be treated as independent bootstrap observations. Coverage and mapping
exclusions accompany every result. No retuning, retry selection, universal
superiority claim or controlled timing claim is authorized by this panel.
