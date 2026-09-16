# Publication Ablation Protocol

Prospective controls for the frozen high-sensitivity/satellite_v2 method.
QfO and OrthoBench remain development-exposed; these experiments explain
components, not independent confirmation. Do not use YGOB outcomes to select
ablation settings or change the frozen validation configuration.

## Factorial Design

Evaluate eight combinations of these binary factors on OrthoBench, then QfO:

- Multi-sequence profile-HMM expansion: the production one-pass expansion
  versus its omission. Retain initial search, RBNH graph, singleton assignment,
  Leiden seed 4, resolution 0.1, and default cluster refinement in both arms.
- Candidate-family expansion: refined seed groups versus satellite_v2
  expanded candidates, with the production expansion parameters unchanged.
- Phylogenetic reconciliation: candidate groups alone versus final root HOGs
  for OrthoBench, and candidate group-derived cross-species pairs versus
  final native ortholog pairs for QfO. The output-level change is part of the
  full pipeline comparison and must be disclosed, not called a pure scorer
  substitution.

In the expansion-plus-phylogeny cells retain satellite_v2's high-confidence
membership constraints. Expansion without reconciliation is a diagnostic of
candidate recall/overmerging, not a proposed deployed mode. Reconstruct merge
traces separately for each seed partition; never transfer constraints from a
different profile-expansion arm.

Infer species trees independently in each end-to-end cell with min_variance
rooting, species_overlap root rule, and positive_paralogy pair rule. A second
diagnostic may hold one frozen inferred tree constant to distinguish tree
estimation effects, but must be labeled a supplied-tree replay.

## Additional Controls

- Preserve pre-refinement and post-refinement partitions for both profile
  arms to measure the separate cluster-refinement contribution. Cluster
  refinement is not synonymous with HMM expansion.
- Add an explicitly labeled unconstrained satellite_v2 reconciliation replay
  to measure the high-confidence membership filter independently.
  Execution specification: use the already-frozen OrthoBench p1_c1 candidate
  partition, replace the membership-constraints argument with the explicit
  --unconstrained-membership diagnostic flag, and retain
  CPU32, the same frozen source/tools and all reconciliation/tree settings.
  Write separate p1_c1_r1_unconstrained_v2 outputs; never overwrite the eight
  factorial cells. Compare to p1_c1_r1 using the existing official score and
  paired reference-family bootstrap (20,000 replicates, seed20260918).
  Treat F1/precision/recall differences as an exploratory three-endpoint
  family with Bonferroni correction. The experiment was named prospectively,
  but this detailed execution specification follows inspection of the core
  factorial results; it is not independent confirmation. It does not use
  YGOB outcomes for selecting a configuration or changing the frozen method.
- Omitting profile expansion does NOT remove the HMM-based initial search.
  A matched sequence-search alternative remains required before attributing
  an advantage to HMMs overall. Freeze its search sensitivity, score
  normalization, candidate coverage, and downstream settings before scoring;
  compare both recovered candidates and runtime rather than assuming equal
  E-values make distinct engines equivalent.

## Reuse And Verification

`replay_high_sensitivity.py` already exports multipass, refined multipass,
profile-expanded, and refined profile-expanded partitions. Its use is gated
on verifying cache/input hashes and reproducing the frozen production
partition. Check gene ordering, self-hits, normalization, significant-hit
filtering, and species-count-dependent refinement before reusing a cache.
For broad panels the production refinement may omit directed-hit arrays;
the generic replay must reproduce that branch before QfO ablations count.

`replay_phylogeny.py` now accepts the production merge-trace checkpoint.
Earlier replays without that input cannot establish a matched satellite_v2
control. They remain exploratory evidence and must not be relabeled.
The restored historical OrthoBench candidate checkpoint validates against
all 8,440 entries of its production merge trace; this validates compatibility,
not a completed matched replay or accuracy result.

Reuse immutable alignment/raw-tree checkpoints only with verified matching
sequences and tool settings. Keep full production runs as timing baselines;
cache replay times are incremental analysis costs, not end-to-end runtimes.
Never overwrite previous output directories or score-based select checkpoints.

## Reporting

Report every cell, neutral/negative effects, coverage, stage times, memory,
and failures. For OrthoBench use the audited weighted statistic, paired
RefOG bootstrap with 20,000 replicates and seed 20260918, and prespecified
factorial contrasts (four conditional effects for each of three factors).
Report nominal intervals and Bonferroni-adjusted intervals across those
12 contrasts times F1/P/R; keep additional diagnostics descriptive.
Recompute the statistic in each replicate, not a mean of per-family F1.
For QfO report all six individual endpoints and coverage; retain its custom
mean as a secondary project summary and do not bootstrap that mean by genes.

Freeze an executable manifest with exact commands and hashes before launching
the cells. This protocol is not evidence that those runs are complete.

## Replay Scheduling Amendment

The label-blind cached OrthoBench equivalence check does not consume YGOB
outputs or require YGOB completion. Its initial exclusive/afterok scheduling
was a resource precaution, not a scientific dependency. To unblock ablation
preparation, run the identical pinned replay command with 32 allocated CPUs
and 64 GiB on a shared node, without the YGOB dependency. Retain the two-hour
limit. Record its time only as an incremental, potentially contended replay
cost, never as controlled end-to-end efficiency evidence. Do not alter YGOB's
allocation or use any held-out scores to select replay settings.

Original job 20919 was confirmed unstarted, then cancelled (zero runtime)
after Slurm refused an in-place sharing update. Its exact submitted command
was recovered with `scontrol write batch_script` and preserved in
`ob_replay_batch_command_20260916.sh`. Scientific source, cache, settings,
inputs, output path and verification criteria are unchanged on resubmission.

## Native Runtime Correction Before Replay V2

Replay 21088 exposed the missing native profile-alignment library; it failed
equivalence and is retained as a defective-runtime diagnostic. Before any new
factorial outcomes, construct a separate checkout `publication_method_native_v2`
of the same full source commit 7f3a9e40dd7e79f842cc2c11fb8b548f9a802806.
Build the three production CPU kernels with GCC 13.3.0 and the Linux flags
used by setup.py (`-O3 -fopenmp -shared -fPIC -march=native`, plus `-mavx2`
for hmm_viterbi). No CUDA binary is present in this CPU-only runtime. This
restores the requested profile stage; no scientific threshold or scoring rule
is changed. Source-only and corrected-native executions are distinct runtimes.

The frozen build manifest `publication_native_runtime_20260916.json` records
source/binary/compiler hashes, complete commands and a passing synthetic
profile probe. The current replay launcher requires that manifest, verifies
the exact native library set before and after inference, and reruns the
exact-interpreter profile probe before creating its output directory.

Run the same historical stage comparison in a new directory
`benchmarks/results/publication_ob_replay_check_v2`, with the same cache,
FASTAs, CPU=32, 64 GiB, two-hour limit, shared-node scheduling and no accuracy
scoring. Pin the amended launcher before submission. The original evidence
and original checkout remain untouched. Success means all four partitions
match, not merely that the profile probe or native process succeeded. Any
remaining discrepancy must be diagnosed before preparing factorial cells.

## Reconciliation Execution Freeze

Corrected replay 21138 passed all four historical partitions byte-for-byte.
Preparation 21161 completed successfully and produced the committed
`orthobench_factorial_prepared_20260916.json` (SHA256
5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382).
Its eight cell definitions are unchanged. Execute the four R=1 cells using
their recorded commands and per-arm constraints; R=0 outputs remain the
already-prepared partitions, not additional inference jobs.

Use at most two concurrent reconciliation tasks, each with 32 CPUs, 64 GiB
and a 24-hour limit on the shared node. Record incremental resource use and
retain failures. The cell wrapper verifies frozen commands, full core source
sets, FASTA/candidate/trace hashes, successful preparation, replay equivalence
and the pinned tool environment before execution, then rechecks dependencies
afterward. It refuses existing output directories and records complete output
inventories. Pipeline exit and file inventory do not replace native validation
or correct output conversion. No official reference benchmark is passed into
the inference process. Freeze/push the wrapper before submission.
