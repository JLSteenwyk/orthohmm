# OrthoHMM HMM Evidence Campaign

## Recommendation

Keep production defaults unchanged. The bounded investigation implemented and
evaluated all four required mechanisms: sequence-weighted family profiles,
alignment-aware recruitment, local joint-tree rescue and selective HMMER
validation. Three standalone configurations pass their original development
screens, but their effects are mixed and their small OrthoBench F1 gains have
paired-family intervals spanning zero. Neither tested combination improves on
its corresponding weighting-only constituent. No generalizable new default or
universal superiority over OrthoFinder's full pipeline is established.

Retain `henikoff`, `henikoff_cluster` and residue-preserving projection as the
[frozen experimental portfolio](../../hmm_campaign_candidate_freeze_20261020.json).
Do not pick a different configuration for each benchmark, revive stopped
branches, retune the failed combinations, or promote the options into defaults.
This is a complete informative development result, not independent biological
confirmation or a claim that every HMM-based approach has reached a plateau.

## Complete-Vector Comparison

Scores are 0 to 1, rounded to six decimals. OB is OrthoBench weighted group F1;
VGNC, SwissTrees and TreeFam-A are QfO pair F1. EC, GO and FAS are conditional
similarity means, not orthology F1. TK is the supplementary Three Kingdoms
reference co-membership F1. Do not average these different statistics.

| Tool / Configuration | OB F1 | VGNC F1 | SwissTrees F1 | TreeFam-A F1 | EC Mean | GO Mean | FAS Mean | TK F1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OrthoHMM canonical (experimental) | 0.741061 | 0.902068 | 0.831348 | 0.615778 | 0.965382 | 0.490236 | 0.763725 | 0.872133 |
| OrthoHMM Henikoff | 0.741965 | 0.900909 | 0.829218 | 0.615728 | 0.967684 | 0.489887 | 0.762825 | 0.871787 |
| OrthoHMM Henikoff, cluster N | 0.741947 | 0.901026 | 0.831348 | 0.613231 | 0.967581 | 0.489826 | 0.762468 | 0.869184 |
| OrthoHMM projection repair | 0.741424 | 0.901388 | 0.833513 | 0.616375 | 0.965835 | 0.489590 | 0.763292 | 0.870227 |
| OrthoHMM high sensitivity (retained) | 0.703590 | 0.665933 | 0.685498 | 0.605008 | 0.932038 | 0.472271 | 0.776585 | 0.826309 |
| OrthoHMM phylogeny (retained) | 0.741061 | 0.901690 | 0.833513 | 0.614864 | 0.965650 | 0.490349 | 0.762993 | 0.872133 |
| OrthoFinder 3.1.5 full | 0.727365 | 0.988546 | 0.848413 | 0.791918 | 0.936130 | 0.469548 | 0.691422 | 0.988582 |
| OrthoFinder 3.1.5 sequence-only checkpoint | 0.587060 | 0.142755 | 0.690515 | 0.734407 | 0.753642 | 0.407043 | 0.561753 | 0.989451 |
| SonicParanoid 2.0.9 | 0.467576 | 0.982794 | 0.798459 | 0.771956 | 0.872761 | 0.454385 | 0.736680 | 0.991276 |
| ProteinOrtho 6.3.6 | 0.450573 | 0.954896 | 0.718111 | 0.643187 | 0.963168 | 0.486336 | 0.813595 | 0.930463 |
| FastOMA 0.3.5 | 0.309069 | 0.950096 | 0.780219 | 0.657902 | 0.897168 | 0.436894 | 0.654366 | 0.861655 |
| OrthoMCL 1.4 | 0.550653 | 0.642027 | 0.766610 | 0.753685 | 0.921326 | 0.463755 | 0.733752 | 0.988929 |

The canonical collection is an experimental matched control, not the shipping
default. The retained phylogeny row is a different historical result collection;
do not substitute it for the current matched control. The other tools, the
sequence-only checkpoint and historical OrthoHMM rows are retained comparator
evidence, not fresh reruns in this campaign. OrthoFinder full uses the official
post-phylogenetic groups, not a later debug N0 partition. A sequence-only
clustering-stage checkpoint is not the same workflow as full orthology inference.

The [full score table](scores.tsv) includes all 29 configuration/control/history
rows, including aliases, four interaction controls, both stopped combinations
and failed executions. These are not 29 distinct tools. All 232 endpoint slots
are retained: 125 scored observations and 107 NA entries. NA means unrun or
unadmitted, never zero. QfO/Three Kingdoms extensions of the stopped combinations
are deliberately unrun under the original gates.

## Precision, Coverage And Populations

| OrthoBench Configuration | F1 | Precision | Recall |
| --- | ---: | ---: | ---: |
| OrthoHMM canonical | 0.741061 | 0.817705 | 0.677553 |
| OrthoHMM Henikoff | 0.741965 | 0.818513 | 0.678509 |
| OrthoHMM Henikoff, cluster N | 0.741947 | 0.818470 | 0.678509 |
| OrthoHMM projection repair | 0.741424 | 0.818590 | 0.677553 |
| OrthoFinder 3.1.5 full | 0.727365 | 0.660651 | 0.809066 |

The OrthoBench comparison is a precision/recall tradeoff, not evidence that
OrthoHMM dominates the full OrthoFinder pipeline. The canonical native universes
are 251378 genes/12 species for OrthoBench, 984137/78 for corrected QfO, and
443217/12 for Three Kingdoms. Native ownership is distinct from benchmark
assignment coverage and reference-gene coverage. Unknown coverage stays NA.

The Three Kingdoms reference covers 2035 genes in 255 families. Canonical
reference-gene coverage is 1.0, precision 1.0 and recall 0.773259; this limited
reference co-membership statistic does not establish genome-wide precision.
TreeFam is supporting evidence and Three Kingdoms is supplementary to the
primary OrthoBench/QfO comparisons.

| QfO Configuration | EC Assessed / Numeric | GO Assessed / Numeric | FAS Eligible / Sampled |
| --- | ---: | ---: | ---: |
| Canonical matched legacy control | 118067 / 118067 | 84021 / 84021 | 5953050 / 176307 |
| Henikoff | 116182 / NA | 84223 / NA | 5950178 / 176709 |
| Henikoff, cluster N | 116167 / NA | 84443 / NA | 5954794 / 175092 |
| Projection repair | 117194 / 117194 | 84269 / 84269 | 5944182 / 174361 |

These conditional means use changing populations; higher means are not
matched-population causal improvements or substitutes for orthology accuracy.
FAS eligible pools are not sample sizes. Saved FAS samples are unseeded and
contain repeated proteins/families. Missing historical numeric counts stay NA.
The complete [endpoint diagnostics](endpoint_diagnostics.tsv) and
[historical observations](historical_cell_diagnostics.tsv) retain precision,
recall, native/reference/assignment coverage, assessed populations, numeric
counts where retained, and resolvable pointers to their actual source records.
Functional similarity endpoints have no invented orthology precision or recall.

## Mechanisms And Stopping Decisions

| Mechanism / Configuration | OrthoBench F1 | Decision |
| --- | ---: | --- |
| Henikoff relative weights, actual N | 0.741965 | Full-vector development pass; mixed effects; experimental |
| Henikoff weights, cluster effective N | 0.741947 | Full-vector development pass; mixed effects; experimental |
| Species-normalized Henikoff, cluster N | 0.737980 | Stop: original precision/recall/F1 gate fails |
| Strict alignment-span guard | 0.739150 | Stop: original F1 gate fails |
| Local joint-tree rescue | 0.741061 | Stop: all 70 reference-family scores unchanged |
| Projection repair alone | 0.741424 | Full-vector development pass; mixed effects; experimental |
| Repaired builtin anchor reselection | 0.736212 | Stop: recall decline exceeds original bound |
| HMMER measurement only, repaired profiles | 0.741424 | Measurement control; identical to repair alone |
| Repaired HMMER anchor reselection | 0.736128 | Stop: recall decline exceeds original bound |
| Henikoff plus projection repair | 0.741784 | Stop: below Henikoff alone |
| Clustered Henikoff plus projection repair | 0.741628 | Stop: below clustered Henikoff alone |

The weighting ablation separates relative sequence weights from effective N;
the cluster count uses the frozen aligned-identity rule, not HMMER's default
entropy weighting. These simplified profiles are not Plan7 models. Every arm
uses its own self thresholds and existing sequence-supported recruitment rules.

Alignment evidence counts known aligned match residues separately from bounding
spans and hit counts. Fragments and asymmetric fusion/fission support are not
rejected merely for different lengths or domain architecture. Local rescue uses
reciprocal profile evidence, constrained joint alignments/trees, reconciliation
and explicit pair/membership checks before ordinary final inference. All 128
contexts and six provisional admissions were evaluated; provisional pair or
membership changes are not accuracy gains.

Selective HMMER used 128 profiles and completely measured 126 of 128 selected
queries; two incomplete decisions were conservatively preserved. It separates
measurement from anchor reselection and retains the matched builtin decision
control. HMM similarity alone does not establish orthology.

All four interaction controls reproduced their original graph, candidate,
constraint, native HOG/pair/confidence and rooted-tree bytes and unrounded
precision/recall/F1. Each combination built 15901 profiles with its actual
declared joint policies. Both improve over canonical but fail the prespecified
stronger-constituent rule. No gate was relaxed or parameter retuned after results.

## Uncertainty And Confirmation

| Standalone Change | OB F1 Difference, Percentage Points | Paired Percentile Interval |
| --- | ---: | --- |
| Henikoff | 0.090379 | [-0.050123, 0.331771] |
| Henikoff, cluster N | 0.088610 | [-0.062645, 0.335866] |
| Projection repair | 0.036341 | [-0.036216, 0.154189] |

The saved procedures use 100000 paired resamples of the 70 reference families.
The original standalone intervals retain their post hoc, development-exposed
qualification. Raw paired and phase-specific Bonferroni-adjusted intervals are
in [the summary](summary.json), including the unchanged
[alignment-guard interval](../aligned_profile_guard_admission_20261018.json).
The interaction preserves six contrasts and
its registered 18-metric adjustment. These are not campaign-wide selection
corrections, proof of equivalence/inferiority, or independent biological
confirmation. Shared family history and fused predictions can violate the
exchangeability assumption. No IID gene-pair or functional-sample uncertainty
claim is substituted for family-level uncertainty.

The three-candidate portfolio was frozen and remotely verified before the
[bounded post-freeze review](../hmm_confirmation_review_20261020.json). Four
tested OrthoBench2 URLs returned 404; two CGOB HTTPS probes failed to connect.
No response bodies or new reference assignments were read. Existing ParaBench
near-total QfO overlap and CGOB sequence-overlap/reference-semantics limitations
preclude a strict independence claim for those staged candidates. This does not
prove that all public benchmarks are unavailable. Missing admitted independent
confirmation limits claims and is nonblocking; no favorable remainder, new
simulator seed or metadata response is treated as biological confirmation.

## Resources, Reproduction And Delivery

The [resource audit](resource_audit.json) identifies the local Threadripper PRO
7995WX and verifies that the 13 full-proteome 16CPU/128G allocations were
sequential, including preserved failed attempts. Other user-account Slurm
allocations were present at every scientific start. These allocation counts
are not instantaneous CPU consumption or an exhaustive host monitor. The DGX
and a quiet window were not required; unrelated analyses were not interrupted.

The [106 phase-cost rows](resources.tsv) preserve wall time, CPU time, measured
maximum-process RSS and separately qualified parent metadata. They are shared-
host cached incremental costs, not fresh end-to-end tool runtimes. Do not sum
RSS, arm totals and their components, or infer a tool speed ranking. Local final
replay excludes earlier nomination/joint-tree preparation. Allocation is not a
measured memory peak. Read-only compilation used 1CPU/16G; the small mechanistic
synthetic diagnostic used 4CPU/16G and is not a full-proteome accuracy benchmark.

The broad regression passed 1647 tests in 43 goal-related modules. The reporting
suite passed 163 tests, including 36 dedicated consolidation cases; these counts
overlap and must not be added as independent coverage. Tests include real CPU
profile scoring, serial/spawn equivalence, exact default arrays, meaningful
384-profile compatibility fixtures, alignment-score equality, actual MAFFT/
FastTree workflows, HMMER 3.4, ownership guards and failure preservation.

The [summary](summary.json) pins the 236636436-byte full source report retained
outside Git. Frozen compiler `fcd95d092cfa72fc109aa51fba637c2af5d401e5` verified
715967 source/content records and all scored diagnostic pointers. The four TSV
archives are byte-identical copies. Its [actual batch script](compiler.batch.sh),
[source compiler](../../report_hmm_campaign_plateau.py),
[protocol](../../hmm_evidence_campaign_20261018.json),
[individual admission](../hmm_individual_effects_20261019/summary.json) and
[interaction admission](../weighted_projection_completed_20261020/summary.json)
retain exact commands, versions, seeds, failure histories and recovery paths.
Earlier failed executions remain distinct from biological gate failures; the
archive-only optional-field correction preserved the same score/cost bytes.
Large proteomes, caches and reports remain outside Git with locations/checksums.

To repeat only the read-only consolidation on this staged workspace, use the
frozen checkout and a fresh output directory:

```sh
env CUDA_VISIBLE_DEVICES= OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONHASHSEED=0 \
  /home/bizon/anaconda3/bin/python \
  benchmarks/work/hmm_campaign_plateau_v1_20261020/benchmark_tools/report_hmm_campaign_plateau.py \
  --repo "$PWD" --output benchmarks/results/hmm_campaign_plateau_reproduction
```

The [completion audit](completion_audit.json) maps the original objective and
prespecified work to authoritative evidence. This delivery closes the bounded
investigation with unchanged defaults and documented limitations; it does not
claim an improved production default, a publication release or independent
superiority over OrthoFinder.
