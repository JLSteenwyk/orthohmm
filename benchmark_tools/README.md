# Benchmark Workflows

## FastOMA comparator

Run FastOMA 0.3.5 on the complete 78-proteome QfO 2020 input and chain the
result into the retained six-metric assessment with:

```bash
benchmark_tools/submit_fastoma_qfo.sh
```

The inference runner uses the cached LUCA OMAmer database, emits bare UniProt
accessions, and forces FastOMA's pairwise ortholog output for the 78-species
dataset. FastOMA requires a supplied species tree, so this comparator uses the
OrthoFinder 3.1.5 tree inferred from the same input proteomes. This tree does
not contain QfO benchmark labels. The runner records the tree checksum and
refuses to overwrite an existing production result.

## Production OrthoHMM performance harness

Use the production CLI benchmark harness for runtime and memory measurements:

```bash
python benchmark_tools/benchmark_production.py \
  bacterial_scaling/subsets/n005 \
  /tmp/orthohmm_n005 \
  /tmp/orthohmm_n005.metrics.json \
  --cpu 32 \
  --accuracy-profile high_sensitivity
```

The JSON records exact commands, input/output checksums, stage wall time,
sampled process-tree RSS, search candidates, significant hits, graph edges,
and orthogroup counts. The bacterial SLURM runner uses this harness rather
than maintaining a second OrthoHMM implementation.

The harness defaults to `--accuracy-profile standard`. Select
`high_sensitivity` explicitly when measuring the bounded k=4 search,
multi-pass graph inference, and strict profile expansion workflow.

## OrthoBench stage diagnostics

Use the offline diagnostic to measure whether reference relationships survive
search hits, graph construction, clustering, and refinement:

```bash
python benchmark_tools/orthobench_stage_diagnostics.py \
  --refogs /path/to/Open_Orthobench/BENCHMARKS/RefOGs \
  --hits-pickle mc100_hits=benchmarks/results/hits_BLOSUM62_mc100.pkl \
  --pairs production_edges=/path/to/orthohmm_edges.txt \
  --clusters production_final=/path/to/orthohmm_orthogroups.txt \
  --official-benchmark /path/to/Open_Orthobench/BENCHMARKS/benchmark.py \
  --json /path/to/stage_diagnostics.json
```

The tool records input hashes and reports both direct reference-pair recovery
and recovery through connected components. Reference labels are read only by
this post-run process; they are never passed into OrthoHMM inference. Pickle
input is provided only for trusted local historical caches.

## High-sensitivity graph replay

Validate production graph and refinement code against a trusted normalized-hit
checkpoint without rerunning the all-to-all sequence search:

```bash
python benchmark_tools/replay_high_sensitivity.py \
  --hits-pickle benchmarks/results/hits_BLOSUM62_mc100.pkl \
  --fasta-directory /path/to/Open_Orthobench/BENCHMARKS/Input \
  --output-directory /tmp/orthohmm_high_replay \
  --official-benchmark /path/to/Open_Orthobench/BENCHMARKS/benchmark.py \
  --json /tmp/orthohmm_high_replay.json
```

The replay uses no reference labels during inference. Its JSON records source
and input checksums, graph counts, runtime, peak process RSS, output checksums,
and optional post-inference OrthoBench scores. Omit `--fasta-directory` to
validate only RBNH, singleton reassignment, and production refinement.

Use `--profile-iterations 2` only to evaluate iterative profile rebuilding.
The production-equivalent default remains one pass. The recorded two-pass
OrthoBench experiment doubled profile-search work and reduced refined F-score
from 70.4 to 70.2, so iterative refinement is not a production recommendation.

Use `python -m benchmark_tools.selective_jackknife_holdout` to test profile
threshold calibration in the presence of mixed, duplicated-species profiles.
The replay flag `--jackknife-single-copy-profiles` applies the jackknife only
when every profile member comes from a distinct species. This mode is separate
from the rejected global `--jackknife-profile-thresholds` ablation.

## Phylogeny-stage replay

Replay the optional phylogeny stage from a frozen candidate-cluster file:

```bash
python benchmark_tools/replay_phylogeny.py \
  --fasta-directory /path/to/proteomes \
  --candidate-clusters /path/to/orthohmm_edges_clustered.txt \
  --output-directory /tmp/orthohmm_phylogeny_replay \
  --json /tmp/orthohmm_phylogeny_replay.json \
  --species-tree-mode infer --cpu 32 \
  --root-rule supported_children --pair-rule positive_paralogy \
  --official-benchmark /path/to/OrthoBench/benchmark.py
```

The JSON records source/input checksums, process-tree memory, stage runtime,
reconciliation counts, output checksums, and optional post-run scores. This is
a development harness; final claims still require fresh end-to-end runs.
For a reconciliation-only comparison, pass `--checkpoint-source` with a prior
replay output. The harness hard-links immutable raw-tree checkpoints when the
filesystem permits, recomputes rooting and reconciliation in a distinct output
directory, and records the source path and checkpoint-hit counts.

## Permissive HMM candidate graphs

Audit a frozen candidate partition before spending time on gene-tree
inference. The diagnostic verifies that candidates only combine immutable seed
families and reports the reference-pair recall ceiling, reference contamination,
family growth, and per-RefOG component recovery:

```bash
python benchmark_tools/candidate_family_diagnostics.py \
  --seed-clusters /path/to/seed_clusters.txt \
  --candidate-clusters /path/to/candidate_clusters.txt \
  --refogs /path/to/OrthoBench/RefOGs \
  --json /tmp/candidate_audit.json
```

Build candidate superfamilies from a trusted normalized-hit checkpoint before
testing phylogenetic reconciliation:

```bash
python benchmark_tools/build_permissive_candidates.py \
  --hits-pickle benchmarks/results/hits_BLOSUM62_mc100.pkl \
  --output-directory /tmp/orthohmm_candidates \
  --json /tmp/orthohmm_candidates/result.json \
  --threshold-factor 0.8 --inflation 1.2 \
  --mcl /path/to/mcl --cpu 32 \
  --official-benchmark /path/to/Open_Orthobench/BENCHMARKS/benchmark.py
```

The support factor scales each gene's weakest reciprocal-best normalized-hit
threshold. Graph construction is label-blind, MCL receives binary edge weights,
and omitted isolates are restored before the candidate file is finalized. The
JSON records parameters, checksums, stage runtimes, peak process-tree memory,
graph/family counts, and any post-inference OrthoBench score.

To preserve OrthoHMM's precise seed groups and merge only mutually
best-supported group pairs, run:

```bash
python benchmark_tools/merge_candidate_superfamilies.py \
  --hits-pickle benchmarks/results/hits_BLOSUM62_mc100.pkl \
  --seed-clusters /path/to/orthohmm_edges_clustered.txt \
  --output-directory /tmp/orthohmm_superfamilies \
  --json /tmp/orthohmm_superfamilies/result.json \
  --max-family-genes 500 --iterations 2 \
  --official-benchmark /path/to/Open_Orthobench/BENCHMARKS/benchmark.py
```

This alternative aggregates directed HMM evidence between seed groups and
joins only reciprocal-best partners. It is intended to avoid giant components
while recovering split remote homologs for subsequent reconciliation. Bounded
iterations let a newly merged group recruit another mutually best-supported
fragment without relaxing the evidence rule.

The experimental `--method satellite` mode relaxes mutual-best ranking only
for a bounded small source group. It still requires bidirectional HMM support,
a source-side best-target margin, size and component caps, and an optional
species-overlap limit. This is intended for fragments whose larger anchor has
another stronger bridge:

```bash
python benchmark_tools/merge_candidate_superfamilies.py \
  --method satellite --hits-pickle /path/to/hits.pkl \
  --seed-clusters /path/to/orthohmm_edges_clustered.txt \
  --output-directory /tmp/orthohmm_satellite_candidates \
  --json /tmp/orthohmm_satellite_candidates.json \
  --max-satellite-genes 5 --max-satellite-ratio 0.75 \
  --min-margin 1.5 --max-satellites-per-anchor 2 \
  --min-coverage 0.5 --min-normalized-support 0.02
```

Validate that bounded behavior independently with:

```bash
python benchmark_tools/candidate_superfamily_holdout.py \
  benchmark_tools/results/candidate_superfamily_holdout.json
```

The synthetic graph contains two-to-five-fragment true families and a cycle of
weaker asymmetric cross-family decoys. It reports pair precision and recall for
one pass and four bounded passes across fixed holdout seeds.

Root-HOG duplication evidence can be replayed without rebuilding alignments or
trees:

```bash
python benchmark_tools/replay_root_hog_rules.py \
  --fasta-directory /path/to/proteomes \
  --candidate-clusters /path/to/candidates.txt \
  --phylogeny-output /path/to/run/orthohmm_phylogeny \
  --output-directory /tmp/root_hog_rules \
  --json /tmp/root_hog_rules/result.json \
  --official-benchmark /path/to/Open_Orthobench/BENCHMARKS/benchmark.py
```

The production default remains `supported_children`; the replay compares that
rule with progressively broader root-level duplication evidence from the same
rooted trees.

Before tuning candidate-family or reconciliation parameters, freeze the
OrthoBench development and validation RefOGs without reading their members:

```bash
python benchmark_tools/create_orthobench_partition.py \
  --refogs /path/to/Open_Orthobench/BENCHMARKS/RefOGs \
  --json benchmark_tools/results/orthobench_partition.json
```

Score a prediction on only one frozen partition with the official
equal-RefOG pairwise formula:

```bash
python benchmark_tools/score_orthobench_partition.py \
  --predictions /path/to/root_hogs.txt \
  --refogs /path/to/Open_Orthobench/BENCHMARKS/RefOGs \
  --partition-json benchmark_tools/results/orthobench_partition.json \
  --partition development --json /tmp/development_score.json
```

Synthetic phylogeny holdouts accept independent reconciliation policies. For
example, evaluate the policy that preserves pairs at mapping-only conflicts
while retaining the conservative root-HOG rule with:

```bash
python benchmark_tools/phylogeny_holdout.py /tmp/phylogeny_holdout.json \
  --work-directory /tmp/phylogeny_holdout \
  --aligner /path/to/mafft --tree-builder /path/to/FastTree --cpu 4 \
  --root-rule supported_children --pair-rule positive_paralogy
```

After predictions are frozen, stratify reconciliation changes by RefOG copy
number and alignment-derived sequence identity:

```bash
python benchmark_tools/analyze_phylogeny_changes.py \
  --refogs /path/to/Open_Orthobench/BENCHMARKS/RefOGs \
  --fasta-directory /path/to/Open_Orthobench/BENCHMARKS/Input \
  --baseline-clusters /path/to/orthohmm_edges_clustered.txt \
  --phylogeny-root-hogs /path/to/orthohmm_root_hogs.tsv \
  --alignment-directory /tmp/refog_alignments \
  --aligner /path/to/mafft --cpu 4 \
  --json /tmp/phylogeny_changes.json
```

This diagnostic reads benchmark labels only after inference. It verifies that
root HOGs remain subsets of their source families, distinguishes splitting
from cross-family merging, counts retained and removed reference relations,
and reports copy-number and median pairwise-identity strata.

## OrthoFinder comparator

This directory pins the OrthoFinder comparator to version 3.1.5 and keeps new
results separate from the historical 2.5.5 outputs.

## Installation

```bash
benchmark_tools/install_orthofinder.sh
```

The installer verifies the upstream archive checksum and the installed version.
Set `SOFTWARE_ROOT` or pass an explicit installation directory when the shared
workspace layout is unavailable.

## Benchmarks

Submit the serialized OrthoBench and QfO workflow with:

```bash
benchmark_tools/submit_orthofinder_v3_benchmarks.sh
```

The QfO runner writes both raw native OrthoFinder pairs (`pairs.tsv`) and pairs
whose IDs occur in the QfO mapping (`pairs.qfo.tsv`). Only the filtered file is
submitted to the official assessment pipeline. The scoring workflow retains
VGNC, SwissTrees, TreeFam-A, EC, GO, and FAS.

To isolate the contribution of OrthoFinder's phylogenetic stages, convert the
pre-phylogenetic DIAMOND/MCL checkpoint retained in `WorkingDirectory`:

```bash
python benchmark_tools/orthofinder_mcl_to_orthogroups.py \
  WorkingDirectory/clusters_OrthoFinder_I1.2.txt_id_pairs.txt \
  WorkingDirectory/SequenceIDs.txt \
  orthogroups_sequence_only.txt
```

This output precedes MSA, gene-tree, species-tree, reconciliation, and HOG
delineation. For QfO it is converted to cross-species co-cluster pairs with
`qfo_benchmark/og_to_pairwise.py`, the same representation used for every
group-only method.

Useful environment overrides include `ORTHOFINDER_BIN`, `SOFTWARE_ROOT`,
`ORTHOBENCH_ROOT`, `QFO_PIPELINE`, `QFO_MAPPING`, `NEXTFLOW_BIN`, and
`JAVA_HOME_17`. Set `RESUME=1` when retrying an interrupted QfO assessment so
Nextflow reuses completed tasks instead of deleting its work directory.

Summarize completed QfO assessments with:

```bash
python benchmark_tools/qfo_summarize_scores.py \
  qfo_benchmark/scoring/orthofinder_v3_diamond
```

Run the production OrthoHMM high-sensitivity QfO workflow separately with:

```bash
benchmark_tools/submit_orthohmm_high_sensitivity_qfo.sh
```

The inference job uses `benchmark_production.py`, records production stage
metrics and manifests, converts final orthogroups to cross-species QfO pairs,
and filters identifiers against the official mapping. Official six-metric QfO
scoring starts only after successful inference.

After freezing a phylogeny-aware candidate, submit its QfO workflow with:

```bash
benchmark_tools/submit_orthohmm_phylogeny_qfo.sh
```

This run freezes `high_sensitivity` plus the `satellite_v2` HMM candidate
profile, minimum-variance internal species-tree rooting, the `species_overlap`
root-HOG rule, and the `positive_paralogy` pair rule. It infers its species tree
internally and submits OrthoHMM's native reconciled cross-species pairs to QfO.
It does not expand root HOGs into cliques.

If a fresh production run reaches a verified candidate-cluster checkpoint but
fails during phylogeny, `submit_resume_orthohmm_phylogeny_qfo.sh` resumes only
that stage. The restart harness verifies and preserves the failed-run
checkpoint, retains both metrics files, and emits `result.json` with combined
wall time and peak process-tree memory before submitting QfO scoring.

Run the secondary Three Kingdoms validation with:

```bash
benchmark_tools/submit_orthohmm_high_sensitivity_three_kingdoms.sh
```

This uses the same production benchmark harness and accuracy profile, then
scores the final orthogroups against the local BUSCO-derived reference. The
workflow expects exactly 12 FASTA files and refuses to replace an existing
result directory.
