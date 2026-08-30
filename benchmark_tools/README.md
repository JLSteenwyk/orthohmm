# Benchmark Workflows

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

Evaluate reciprocal profile-HMM repair without final benchmark labels with:

```bash
python -m benchmark_tools.reciprocal_profile_merge_holdout \
  --development-seed 1103 \
  --holdout-seed 2909 \
  --families 36 \
  --cpu 8 \
  --output /tmp/reciprocal_profile_merge_holdout.json
```

The harness splits each synthetic family into two three-species profiles and
includes close paralogs as negative merge candidates. The replay flag
`--reciprocal-profile-merges` enables the independently screened rule; it is
off by default until fresh production benchmarks pass.

Use `python -m benchmark_tools.selective_jackknife_holdout` to test profile
threshold calibration in the presence of mixed, duplicated-species profiles.
The replay flag `--jackknife-single-copy-profiles` applies the jackknife only
when every profile member comes from a distinct species. This mode is separate
from the rejected global `--jackknife-profile-thresholds` ablation.

The same harness tests `--profile-score-per-target-residue`, which compares
local HMM scores and their member-derived thresholds after dividing each by the
target protein length. Use both screened flags together for the fragment-aware
candidate:

```bash
python -m benchmark_tools.selective_jackknife_holdout \
  --development-seed 1103 \
  --holdout-seed 2909 \
  --families 36 \
  --cpu 4 \
  --output /tmp/profile_length_calibration_holdout.json
```

On both synthetic splits, per-target-residue calibration retains perfect
close-paralog and unrelated rejection while recovering partial proteins missed
by raw full-length score thresholds. Both replay flags remain off by default
until fresh production benchmarks pass.

Evaluate direct profile-profile repair with:

```bash
python -m benchmark_tools.profile_profile_holdout \
  --development-seed 1103 \
  --holdout-seed 2909 \
  --families 36 \
  --cpu 4 \
  --output /tmp/profile_profile_holdout.json
```

`--profile-profile-merges` aligns only exact-4-mer-prefiltered profile pairs,
normalizes each local alignment by both profile self-scores, and connects only
mutual nearest profiles above 0.60. Candidate capacity is 30 per profile and a
pair may contain at most 80 seed genes, bounding both alignment work and graph
growth. The mode is disabled by default: 0.60 was inactive on OrthoBench, while
the independently screened 0.40 sensitivity follow-up reduced F-score.

Evaluate HMM-supported component splitting with:

```bash
python -m benchmark_tools.profile_component_split_holdout \
  --development-seed 1103 \
  --holdout-seed 2909 \
  --families 36 \
  --cpu 4 \
  --output /tmp/profile_component_split_holdout.json
```

`--component-split-high-duplication` retains separate strong-edge components
inside paralog-rich clusters before falling back to the historical single-core
degree rule. Profile-HMM anchor edges participate in these components, while
edges below the existing 1.5 component threshold do not. The mode is disabled
by default pending fresh production benchmarks.

Evaluate direct HMM fallback for pairwise-unanchored candidates with:

```bash
python -m benchmark_tools.direct_profile_fallback_holdout \
  --development-seed 1103 \
  --holdout-seed 2909 \
  --families 36 \
  --cpu 4 \
  --output /tmp/direct_profile_fallback_holdout.json
```

`--direct-profile-fallback` requires single-copy jackknife calibration. It adds
one confidence-ratio edge only when a significant calibrated profile winner has
no positive pairwise hit to any member of its selected cluster. Edge weights
are bounded to the existing graph scale, and the mode is disabled by default.

Evaluate calibrated pair-seeded profiles with:

```bash
python -m benchmark_tools.pair_profile_holdout \
  --development-seed 1103 \
  --holdout-seed 2909 \
  --families 36 \
  --cpu 4 \
  --output /tmp/pair_profile_holdout.json
```

Set `--profile-min-cluster-size 2 --pair-profile-threshold-ratio 0.7` together
with `--jackknife-single-copy-profiles` to include cross-species two-gene seed
clusters. Their thresholds come from both leave-one-out directions using a
duplicated retained-member profile; same-species pairs remain at the strict
in-sample threshold. Pair profiles remain disabled by default pending a fresh
production benchmark.

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

Run the secondary Three Kingdoms validation with:

```bash
benchmark_tools/submit_orthohmm_high_sensitivity_three_kingdoms.sh
```

This uses the same production benchmark harness and accuracy profile, then
scores the final orthogroups against the local BUSCO-derived reference. The
workflow expects exactly 12 FASTA files and refuses to replace an existing
result directory.
