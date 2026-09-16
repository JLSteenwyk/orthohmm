# YGOB Group-Recovery Validation Protocol

Scientific specification frozen before inference or accuracy inspection.
Input preparation and exact-sequence overlap counts were inspected; no method
predictions on this dataset have been examined. This is a novel-taxon,
curated-reference transfer experiment, not a family-disjoint claim.

## Inputs And Reference

- Use the acquired YGOB v7 snapshot identified by `ygob_overlap_20260916.json`.
- Exclude the complete Saccharomyces genus: Scerevisiae, Skudriavzevii,
  Smikatae, Suvarum. Retain the other 16 named species in
  `ygob_validation_inputs_20260916.json`; do not substitute species after scoring.
- Exclude OFF proteins. Uppercase sequences and remove one terminal stop.
  Exclude internal-stop proteins; after genus exclusion none remain.
  Reject any other unsupported residue rather than silently changing sequences.
- Retain 83,404 proteins in inference. The prepared FASTAs preserve original
  IDs without using reference labels in the inference inputs.
- Exclude both reference rows containing duplicate membership (113 and 9896).
  Their 13 retained proteins remain in inference but are outside the scoring
  universe. There are 83,391 scored genes in 10,250 pillars, including singletons.
- The target is recovery of curated homolog groups spanning pre- and post-WGD
  species. It is not resolved pairwise orthology. YGOB's two post-WGD columns
  are not A/B assignments; neither column order nor pillar cliques may be
  used as ground truth for a resolved-ortholog claim.

## Frozen Methods

OrthoHMM production code is pinned to `7f3a9e4` in a detached worktree.
Relative to the historical satellite_v2 inference commit, the only core
change is preserving a candidate checkpoint; no new inference policy is
introduced. The existing `high_sensitivity` profile includes its fixed Leiden
seed of 4. `PYTHONHASHSEED=0`; BLAS/OpenMP environment threads are limited to one.

1. OrthoHMM high_sensitivity: built-in HMM-centered pipeline, Leiden,
   BLOSUM62, E-value 1e-4, CPM resolution 0.1, default refinement profile,
   CPU 32 and eight threads per worker; no phylogeny.
2. OrthoHMM satellite_v2: the same base settings, candidate profile satellite_v2,
   internally inferred species tree, min_variance rooting, species_overlap
   root duplication rule, positive_paralogy pair rule. Evaluate final root
   HOGs, not its pairwise orthology file, against this group reference.
3. OrthoFinder 3.1.5 full: DIAMOND, `-t 32 -a 8 -S diamond`, other settings
   unchanged from its installed version. Evaluate finalized orthogroups.
4. OrthoFinder sequence-only: extract the MCL checkpoint from that same run
   using the tested converter. This is a diagnostic checkpoint, not an
   independently timed full run or a second root-HOG prediction.

Exact inference commands are in `run_ygob_validation.slurm`. That workflow
checks input hashes, refuses to overwrite existing outputs, records source
and environment provenance, and performs no accuracy scoring. Request an
exclusive scheduler allocation and run methods sequentially. Verify actual
node workloads before treating measurements as controlled efficiency results;
exclusive scheduling alone does not exclude unscheduled processes. GNU-time
RSS and OrthoHMM process-tree RSS must remain separately labeled.

## Endpoints And Analysis

- Primary contrast: satellite_v2 versus full OrthoFinder in micro-averaged
  unordered gene-pair F1 for group co-membership within the 83,391-gene scoring
  universe. Include within-species pairs: the target groups predate the WGD.
  Every predicted gene outside this universe is removed for this statistic.
- Report precision, recall, reference-gene coverage, predicted-group coverage,
  and exact-group recovery alongside F1. Singleton reference groups contribute
  no true pairs but their erroneous merging contributes false-positive pairs.
- Secondary contrast: high_sensitivity versus full OrthoFinder. The sequence-only
  checkpoint is a diagnostic of phylogenetic processing, not a primary contrast.
- Paired reference-pillar bootstrap, 20,000 replicates, PCG64 seed 20260917.
  Recompute F1/P/R from sufficient counts in each replicate. Allocate half of
  each cross-reference false-positive pair to each incident reference pillar
  so their summed contribution equals the actual micro statistic.
- Report nominal 95% percentile intervals and Bonferroni-adjusted intervals
  across six contrasts/metrics (two OrthoHMM modes times F1/P/R). These are
  approximate and assume exchangeable pillars; no gene-pair independence.
- Retain all results, including negative effects and failed inference. Do not
  alter inputs, parameters, endpoints, or exclusions after accuracy inspection.
  A method change based on these outcomes requires new independent confirmation.

## Overlap Gate Before Accuracy Inspection

Complete a label-independent homology screen against the development inputs
before examining accuracy results. The screen is descriptive, not a means
to select advantageous families. Use DIAMOND 2.1.11 very-sensitive search,
E-value <=1e-5, identity >=30%, query and subject coverage >=50%, retaining
the best qualifying hit per query. Record exact command, input hashes, and
fractions of candidate proteins/pillars with qualifying hits. Do not infer
absence of homology from a negative screen.

Check shared reference resources separately. Shared gene families with QfO
are expected, and are a limitation even when candidate taxa and curated
labels have not been used for method selection. No result from this protocol
alone can establish unrestricted generalization to arbitrary datasets.

## Completion Gates

Before interpreting this experiment: verify the overlap report, scored
reference construction, matched input manifests, tool completion and versions,
native-output conversions, independent score arithmetic, and paired uncertainty.
Further biological case studies, HMM ablations, simulation robustness,
reference-resource auditing, and manuscript work remain required by the full
publication goal. This protocol does not replace those requirements.

## Native Runtime And Scheduling Amendment

Original job 20917 was cancelled while still pending (zero runtime), after
label-blind OrthoBench replay revealed that its source-only checkout lacked
the profile-alignment binary. No YGOB outcome has been inspected. Resubmit
using separate CPU-native checkout `publication_method_native_v2`, retaining
the exact source commit, prepared inputs/reference, scientific settings and
32-CPU/eight-worker-thread allocation. Corrected cached replay 21138 now
matches all four historical OrthoBench partitions byte-for-byte.

The corrected launcher requires a pinned launcher revision and build-manifest
hash, verifies all native binaries, and records exact-checkout profile probes
before and after both OrthoHMM stages. It preserves copies of the submitted
launcher and runtime manifest in the result directory so provenance does not
depend on a disappearing Slurm spool file. Pin MAFFT/FastTree/DIAMOND lookup
to the benchmark installations, record entrypoint hashes and both Python
package inventories, and verify OrthoFinder 3.1.5 before inference.

Use a shared-node allocation (32 CPUs, 128 GiB, 24-hour limit), rather than the
original exclusive request, to allow accuracy validation to progress alongside
other work. This is a prospective scheduling deviation, not controlled timing
evidence. Retain per-method GNU-time logs but label them contended; matched
efficiency measurements remain a separate unfulfilled requirement. Do not
change scientific settings or use test outcomes to select a runtime.
