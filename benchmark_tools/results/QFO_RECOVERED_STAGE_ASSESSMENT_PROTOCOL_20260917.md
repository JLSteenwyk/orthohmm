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

## Frozen Assessment Inputs

Conversion job21522 completed0:0 in2:45. The
[four-stage preparation report](qfo_recovered_stage_pairs_20260917.json) has SHA256
ce6f19cd005b886a91dc3aee63cb7b41e7f448d8cad5e65ff7cb10d92a636100.

| Stage | Raw pairs | Retained pairs | Mapping exclusions |
| --- | ---: | ---: | ---: |
| Multipass | 32102853 | 32012674 | 90179 |
| Multipass refined | 8741224 | 8710340 | 30884 |
| Strict profiles | 31768877 | 31682321 | 86556 |
| Strict profiles refined | 8741773 | 8710722 | 31051 |

These counts are preparation evidence, not accuracy. The
[scoring environment](qfo_assessment_environment_20260917.json) has SHA256
e86545fd04cb644ed642cf2a31b4993fb2225ca4c05743cd03e661db129e82bc.
It records124tracked pipeline files atc0854a96c1a0fd7f2a891d971af0863002fabc90,
103reference files,466Java runtime files, three local2022.1container images,
Nextflow22.10.8build5860, SingularityCE4.3.2 and local configuration/support
executables. This is a present-day identity snapshot, not retrospective proof
of every historical dependency or a hermetic host image.

`qfo_recovered_assessment.config` selects the checksum-pinned local container
paths and a local eight-CPU64GB executor budget. The actual Nextflow effective
configuration was evaluated offline and retained. The new runner checks hashes
before/after each command, refuses existing result/work/launch directories and
never invokes deletion or automatic resume. A four-task array runs at most one
stage at a time with eightCPUs64GiB24hours per stage. Scoring exit0 remains
pending independent validation of allsix native endpoints and their semantics.
