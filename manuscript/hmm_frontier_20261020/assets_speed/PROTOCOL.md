# Small-dataset OrthoHMM generation comparison

This descriptive experiment compares the unchanged original `v0.1.0` tag
(`7d1540c624ab19a9d6083d9bcbdfe55c7fcd26c2`) with a frozen development
implementation named OrthoHMM2 in the manuscript. OrthoHMM1 is a generation
label, not a claim that the old package was released as version 1.0. The exact
modern revision, Python, native kernels, external tools, inputs and observer
are hashed in `plan.json` before timing begins.

The dataset consists of the three smallest complete proteome files, selected
by protein count from the existing bacterial `n005` set before observing
timings. Selection does not use orthology labels or benchmark scores. All
proteins and sequences are retained without downsampling or identifier edits.

Each of three configurations receives three fresh native executions in a
balanced order, serially, in unique output directories. The configurations are
original phmmer/MCL, development high-sensitivity sequence-based inference,
and the latter with inferred phylogenetic refinement. All stop after completed
orthogroup inference, before optional sequence-file exports. The phylogenetic
configuration also produces native trees and ortholog relations, so it performs
additional work. This is not an equal-output or equal-accuracy comparison.

The fixed eligible CPU set is 64-79, sixteen distinct physical cores on the
local Threadripper. The original CLI uses `--cpu 4`, which starts four concurrent
phmmer processes each with four compute threads. Its later graph/MCL phases
retain four workers. The development CLI uses `--cpu 16` with four threads per
worker. Affinity is not an exclusive allocation; other jobs, including jobs on
SMT siblings, can compete. All timings are shared-host observations, not an
isolated machine-speed claim. No unrelated jobs are stopped or suspended.

External monotonic elapsed time includes Python startup and imports. Initial
homology-search times are recorded separately using existing development
metrics and a top-level-call observer around the unchanged original code.
The observer changes no original source or search command. Native process
execution ends before output validation and provenance hashing. OS caches
are not flushed. A 1200-second per-attempt timeout is explicit in the plan;
failed or timed-out attempts cannot become successful runtimes or zeroes.

Admission requires a successful native exit, complete metrics, exact
once-only ownership of all input proteins, and (for the original) all nine
directional phmmer tables ending in the completion marker. Summaries require
all nine planned attempts to pass, retain all observations, and use medians
and observed ranges rather than fastest repeats or confidence intervals.
Original/modern ratios are ratios of elapsed-time medians, descriptive only.

From a checkout carrying the harness, with the recorded scientific source
revision and tools available:

```sh
python -m benchmark_tools.compare_orthohmm_generations prepare --output NEW_RUN \
  --modern-revision RECORDED_MODERN_COMMIT --mcl /path/to/compatible/mcl
python -m benchmark_tools.compare_orthohmm_generations run --output NEW_RUN
python -m benchmark_tools.compare_orthohmm_generations summarize --output NEW_RUN
```

Preparation and execution refuse existing attempt directories. Retain failed
attempts and use a new plan directory for any protocol change. No score or
accuracy-frontier coordinates are updated by this experiment.

The first 22 October pilot found that the shell's default MCL 02-063 rejects
the original CLI's `--abc` option. That failed original run is not a completed
runtime. The pilot controller was stopped before continuing with the invalid
dependency. Its completed and interrupted attempts remain separately recorded.
A compatible MCL 22-282 build was checked using the pilot search tables only
for an explicitly untimed dependency smoke test, then the entire nine-attempt
panel was restarted with new searches and new output directories. The revised
preflight rejects MCL builds without `--abc` support.
