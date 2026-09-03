# OrthoHMM phylogeny hill-climb, 2026-09-03

## Decision

The strict combined target has not been met. The strongest OrthoBench method,
depth-two topology resolution v7, reaches F=74.297 and beats full OrthoFinder
3.1.5 by 1.597 points. The strongest QfO result, sparse-overlap v4, reaches a
six-metric mean of 0.752162, which remains 0.029909 below OrthoFinder's
0.782071. The v4 and v7 QfO results are complete checkpoint previews rather
than fresh end-to-end runs, so neither can count as a final win even if its
score had crossed the target. A fresh production run using the frozen
satellite-v2 settings reaches only 0.748243, confirming the deficit
independently.

Phylogeny should remain opt-in. The current default nonphylogenetic behavior
has not changed. The fresh positive-paralogy run completed; the fresh
supported-paralogy run failed in an isolated Leiden worker before candidate
assembly and was not scored.

## OrthoBench

| Method | Run type | F | Precision | Recall | Exact |
| --- | --- | ---: | ---: | ---: | ---: |
| OrthoFinder 3.1.5 sequence-only MCL | Full-run checkpoint | 58.700 | 45.700 | 82.000 | 18 |
| OrthoHMM high sensitivity | Fresh end to end | 70.400 | 78.900 | 63.500 | 13 |
| OrthoHMM original inferred phylogeny | Fresh end to end | 70.900 | 80.400 | 63.500 | 14 |
| OrthoFinder 3.1.5 full | Fresh end to end | 72.700 | 66.100 | 80.900 | 19 |
| OrthoHMM satellite v2 / sparse v4 | Fresh root HOGs / replay pairs | 74.106 | 81.770 | 67.755 | 15 |
| OrthoHMM depth-two v7 | Frozen replay | **74.297** | **81.565** | 68.218 | 15 |

The v7 split is development F=73.051, validation F=75.332. Method selection
used only the predefined development partition, then opened validation once.
No QfO result was used to select v7 or any v9 profile-expansion candidate.

## QfO

| Method | VGNC F | SwissTrees F | TreeFam-A F | EC | GO | FAS | Mean |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OrthoFinder 3.1.5 sequence-only | 0.14029994 | 0.69950736 | 0.71101320 | 0.75339446 | 0.40760020 | 0.55794138 | 0.54495942 |
| OrthoHMM high sensitivity | 0.66744809 | 0.67385086 | 0.57586996 | 0.93104489 | 0.47228231 | 0.77479433 | 0.68254841 |
| OrthoHMM original inferred phylogeny | 0.74536115 | 0.64514895 | 0.34793251 | 0.97428495 | 0.48504497 | 0.77917375 | 0.66282438 |
| OrthoHMM satellite v2 preview | 0.90053270 | 0.79675097 | 0.57039560 | 0.97337863 | 0.48877334 | 0.76094581 | 0.74846284 |
| OrthoHMM satellite v2 fresh | 0.90073006 | 0.79197597 | 0.57549491 | 0.96919541 | 0.49012048 | 0.76193894 | 0.74824263 |
| OrthoHMM supported-paralogy v3 preview | 0.90463576 | 0.80189627 | 0.57824448 | 0.97260791 | 0.48782482 | 0.76000571 | 0.75086916 |
| OrthoHMM sparse-overlap v4 preview | 0.90645819 | 0.80609201 | 0.58610892 | 0.97083545 | 0.48560258 | 0.75787621 | **0.75216223** |
| OrthoHMM depth-two v7 preview | 0.90651558 | 0.80425740 | 0.58634563 | 0.97074825 | 0.48549269 | 0.75655274 | 0.75165205 |
| OrthoHMM root closure v8 | 0.66149862 | 0.67436355 | 0.57355900 | 0.95679382 | 0.47202445 | 0.75525102 | 0.68224841 |
| OrthoFinder 3.1.5 full | **0.98822936** | **0.85916044** | **0.74288699** | 0.94198234 | 0.46816557 | 0.69200063 | **0.78207089** |

OrthoHMM v4 exceeds OrthoFinder on EC, GO, and FAS, but the deficits in all
three reference-tree tests dominate the mean. It is particularly short on
TreeFam-A sensitivity: TPR=0.426331 versus OrthoFinder's 0.620585, despite a
higher PPV of 0.937437 versus 0.925226.

The authoritative fresh satellite-v2 run predicts 5,483,453 mapped pairs and
has VGNC TPR/PPV=0.819629/0.999643, SwissTrees=0.679322/0.949421, and
TreeFam-A=0.410899/0.960079. Its 0.748243 mean is 0.033828 below OrthoFinder.
Inference used 32 CPUs, 75,848 seconds of pipeline wall time, and 16.752 GiB
peak process-tree RSS. The QfO scorer used 8 CPUs and completed in 1,609
seconds; its largest task was FAS at 34.2 GB peak RSS.

## Why OrthoFinder Leads

OrthoFinder's advantage is not simply that it uses a tree. Its sequence-only
MCL checkpoint is very sensitive but severely over-grouped: it predicts about
163.0 million mapped QfO pairs and scores only 0.544959. The full pipeline
reduces that to 13.4 million pairs while retaining high topology recall and
raising the mean to 0.782071. It starts from a much broader candidate graph,
then uses gene trees and species-tree reconciliation to recover precision.

OrthoHMM starts from cleaner but narrower candidates. Satellite v2 raises
OrthoBench candidate reference-pair recall from 0.495839 to 0.600717 at
0.997229 precision. Sparse pair extraction then gives very high QfO PPV, but
it cannot infer relationships between homologs that never entered the same
candidate family. Once that ceiling is reached, relaxing reconciliation only
creates false pairs.

The topology components make the ceiling explicit:

| Method | VGNC TPR / PPV | SwissTrees TPR / PPV | TreeFam-A TPR / PPV |
| --- | ---: | ---: | ---: |
| OrthoHMM high sensitivity | 0.833 / 0.557 | 0.706 / 0.645 | 0.447 / 0.808 |
| OrthoHMM sparse v4 | 0.836 / 0.990 | 0.704 / 0.943 | 0.426 / 0.937 |
| OrthoHMM root closure v8 | 0.838 / 0.546 | 0.708 / 0.644 | 0.445 / 0.806 |
| OrthoFinder 3.1.5 full | 0.982 / 0.994 | 0.788 / 0.944 | 0.621 / 0.925 |

Root-HOG closure restores roughly the nonphylogenetic TPR ceiling but destroys
the precision gained by reconciliation. This rules out pair closure as the
primary solution and locates the remaining deficit upstream.

## Profile-HMM Iteration

The v9 work introduced an independently configurable, label-blind harness for
rebuilding profile HMMs from seed groups, searching for candidate satellites,
and requiring reciprocal profile coverage. It is retained as research tooling
but is not enabled in production defaults.

| Candidate experiment | Development F | Precision | Recall | Decision |
| --- | ---: | ---: | ---: | --- |
| High-sensitivity seeds | 64.754 | 69.501 | 60.613 | Baseline |
| Rebuilt profile HMMs | 64.921 | 68.770 | 61.480 | Negligible gain |
| Rebuilt profiles then v7 | 70.778 | 84.404 | 60.939 | Below v7 |
| Second pass after satellite v2 | 58.325 | 50.227 | 69.537 | Contaminated |
| Second pass, overlap <=0.5, then v7 | 73.002 | 80.037 | 67.103 | 0.049 below v7 |
| Four-pass reciprocal sequence links | 49.190 | 38.349 | 68.577 | Transitive over-grouping |
| Reduced-alphabet profile prefilter | 64.518 | 68.906 | 60.654 | Slower and below baseline |
| Source groups up to 500 genes | 58.043 | 54.971 | 61.480 | Precision collapse |
| Reciprocal profile coverage >=0.5 | 64.754 | 69.501 | 60.613 | Inactive |

Relaxing the profile score ratio, building profiles from two-member groups,
and unioning normal and reduced-alphabet links produced no score improvement.
Preserving every original seed boundary after root-HOG coalescence regressed to
F=66.581. Unioning all cross-species seed pairs raised development direct-pair
recall from 0.570 to 0.710 without cross-RefOG pairs, but QfO high-sensitivity
already demonstrates that this pair ceiling remains well below OrthoFinder.

## Prefilter Tests

A targeted *C. elegans* to human diagnostic tested the hardest apparent
prefilter hypotheses using development labels only after candidate scoring.
Doubling the k=4 candidate cap from 100 to 200 increased candidates from
2,000,717 to 3,964,268 and retained HMM hits from 44,878 to 52,862. Of the
7,984 additions, only one was within a development RefOG, 33 linked a
reference gene to a nonreference gene, and 7,950 were other links.

Adding exact k=3 candidates contributed 5,482 retained hits but no same-RefOG
or cross-RefOG development link. Relaxing the HMM E-value from 1e-4 through
1e-1 increased retained hits from 44,878 to 66,968 while the number of
same-RefOG links stayed exactly 10. These results reject simple cap, k-mer,
and E-value relaxation as routes to the missing remote homologs.

## Limiting Stage

The evidence establishes a technical blocker for further threshold
hill-climbing within the current exact-k-mer-prefiltered sequence-to-profile
architecture. Candidate generation and family assembly do not expose enough
separable remote-homolog evidence. Reconciliation variants can choose among
existing candidates but cannot recover absent relationships; broad closure or
weak transitive links recover recall only by collapsing precision.

The next justified architecture is a broad weighted candidate source from a
sensitive engine such as MMseqs2 or DIAMOND, or true profile-profile search for
difficult families, followed by profile-HMM membership validation, domain and
length compatibility checks, bounded family assembly, and supported
gene-tree/species-tree reconciliation. HMM evidence should remain central,
while the broad search should be treated only as candidate generation. This
must be selected on synthetic holdouts and OrthoBench development before one
validation and QfO evaluation.

## Production Reliability

Fresh v2 inference job 20873 completed. Its first QfO scoring job, 20874,
failed because Darwin truncated an overlong `predictions.db` path. Job 20883
reran from scratch with `/tmp/qw/work-v2f` and completed all six assessments.

Fresh supported-paralogy job 20880 completed the 18.5-hour search and profile
expansion, then an isolated Leiden process received SIGSEGV. The remaining
cluster file was from an earlier Leiden pass, so it was not accepted as a
candidate checkpoint and dependent scoring job 20881 was cancelled. Commit
`a7134bd` now retries one SIGSEGV, removes stale stage output before each
attempt, and atomically publishes completed clusters. This reliability fix
passed the 394-test unit suite but was made after the run and is not presented
as an accuracy result.

The full source manifests of both long runs match the commits captured at
launch (`681ce931` for v2 and `b4ca0aab` for v3), despite the harness recording
a later repository HEAD and dirty generated outputs when each run finalized.

## Reproducibility

The full v9 ledger, exact commits, run counts, resource use, checksums, and
negative decisions are in
`benchmark_tools/results/phylogeny_iterative_hmm_v9_negative_20260903.json`.
The strongest frozen artifacts are
`phylogeny_depth_two_v7_freeze_20260902.json` and
`qfo_sparse_overlap_v4_checkpoint_preview_20260902.json`. Fresh v2 provenance
and scores are in `qfo_satellite_v2_fresh_20260903.json`; the unsuccessful v3
production attempt is in `qfo_satellite_v3_failed_20260903.json`. No benchmark
taxon, RefOG identifier, or QfO label is present in production inference logic.

Verification: 394 unit tests pass; focused Ruff and scoped whitespace checks
pass for every touched source, test, report, and result file.
The integration suite has eight failures before v9 assertions because the
installed MCL 2.0 binary rejects the modern `--abc` option; generated
`tests/samples/**` outputs were left untouched and uncommitted.
