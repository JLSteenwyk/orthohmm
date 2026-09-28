# Native FAS Annotation Path Panel

This post hoc mechanism diagnostic follows the
[native omission probe](QFO_FAS_OMISSION_MECHANISM_20260928.md). It counts
feature paths for every protein record in all 78 annotation JSON files pinned
by the original QfO assessment environment. No proteomes or proteins are
selected according to method scores. The diagnostic does not rerun pairwise
FAS scoring or modify any original benchmark endpoint.

## Native Computation

The retained image contains greedyFAS 1.18.7. The actual `runMultiTaxa` argument
parser and `calcFAS.fas` option builder supply effective settings, intercepted
only at entry to `fc_start` before pairwise scoring. They preserve the native
feature-type configuration, e-value cutoffs, overlap settings and `10**15`
path limit. The script then calls native feature linearization, region mapping
and path counting without mocks or a replacement graph algorithm.

Each annotation file is checked against the frozen input hash. Per-protein
counts, instance counts and strict-greater-than-limit flags are retained in
local work artifacts; their hashes and summaries appear in the panel report.
The separate completeness reviewer checks that each annotation protein has
exactly one result, validates integer counts and exclusion flags, recomputes
all summaries and checks for identifiers shared across annotation files.
This checks completeness and consistency, not an independent path-counting
algorithm.

The initial yeast pilot counted 6,049 protein records, none above the limit;
its maximum was 56,159,488 paths. All 35 focused record-coverage, panel-binding
and attrition tests pass.

## Complete Results

The [full panel](fas_native_path_panel_20260928.json) completed successfully:
all 78 annotation files contain 984,137 protein records, of which 1,143 exceed
the native path limit (0.1161%). Exclusions occur in 34 annotation files. The
separate [completeness review](fas_native_path_panel_review_20260928.json)
reproduces the totals, verifies every annotation identifier occurs exactly
once in its result, and finds no protein identifier shared between annotation
files. Thus these record totals also equal unique identifier totals here.

The [execution receipt](fas_native_path_panel_execution_20260928.json)
records successful native and review exits, the exact invocation and the
container pin recheck. Per-file hashes bind the complete retained work
artifacts. This is not a second full native computation or an independent
validation of the native graph algorithm.

## Scope

A protein above this limit can trigger the previously verified omission path
when it occurs in a newly scored pair. Its presence does not identify the
historical requested sample or prove why any particular saved score is absent.
Precomputed scores are a separate stratum and must not be relabeled missing
merely because current annotation complexity exceeds the new-scoring limit.
The panel alone estimates neither missing-pair frequency for a method nor a
comparison bias, population mean or confidence interval. Linking method
predictions, precomputed membership and historical sampling remains separate
work.

The workflow uses two worker processes on the shared Threadripper. It is not
part of the controlled 27-run timing panel and supplies no speedup claim.
Large per-protein work artifacts are not committed.

```sh
singularity exec --cleanenv \
  qfo_benchmark/scoring/container_cache/qfobenchmark-fas_benchmark-2022.1.img \
  python benchmark_tools/run_native_fas_path_panel.py \
  --environment benchmark_tools/results/qfo_assessment_environment_20260917.json \
  --directory NEW_WORK_DIRECTORY --output NEW_PANEL.json
python -B -m benchmark_tools.review_fas_path_panel \
  --panel NEW_PANEL.json --output NEW_REVIEW.json
```
