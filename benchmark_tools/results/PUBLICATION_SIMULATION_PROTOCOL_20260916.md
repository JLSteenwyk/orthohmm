# Prospective Simulation Robustness Protocol

Scientific specification recorded before method inference or accuracy
inspection on this panel. The two previous four-species smoke histories
tested infrastructure and truth conversion only; they are not panel members.
This protocol does not establish that panel generation or inference completed.

## Engine And Base Design

Use Zombi `8db13ee4ba007f46c17f38586d31e5aa617c1647`, Pyvolve 1.1.0,
and the recorded 128-bit family-seed adapter from `run_zombi_seeded.py`.
Use the versions validated in `ZOMBI_SMOKE_20260916.md`, pinning all resolved
dependencies and workflow source before launch. Do not tune simulation
parameters to make either method win. Do not use YGOB outcomes to change
this panel or the frozen inference method.

Use ten independent master seeds: 20261001 through 20261010 inclusive.
Each base simulation has eight extant species, 100 root-origin genes,
protein length 300, WAG substitution, speciation rate 1, extinction rate 0,
and lineage-count stopping at eight species. Retain native unscaled time
trees (`SCALE_TREE=0`); protein substitution scaling is a separate parameter.
Set genome minimum size to one, event extensions to geometric parameter one,
and transfer, origination, inversion, and transposition rates to zero.
Request native events, genomes, gene trees, and reconciled XML. All other
native defaults must be expanded into explicit parameter files and hashed
in the executable launch manifest. Unexpected unsupported output or event
types must fail validation rather than be silently excluded.

## Conditions

| Condition | Genome duplication rate | Genome loss rate | Sequence scaling | Post-simulation operation |
| --- | ---: | ---: | ---: | --- |
| baseline | 2 | 1 | 0.2 | None |
| divergent | 2 | 1 | 0.8 | None |
| turnover | 10 | 8 | 0.2 | None |
| divergent_turnover | 10 | 8 | 0.8 | None |
| missing20 | Reuse baseline | Reuse baseline | Reuse baseline | Remove 20% of genes independently of labels |
| uneven_taxa | Reuse baseline | Reuse baseline | Reuse baseline | Thin a tree-defined clade |
| taxon_count_control | Reuse baseline | Reuse baseline | Reuse baseline | Hash-sampled taxa, count matched to uneven_taxa |

Rates above are genome-level event rates in Zombi G mode, not per-gene
rates. Record realized duplications/losses, surviving copies, total branch
length, and observed sequence-divergence distributions. Labels such as
"turnover" name designed conditions, not a guarantee that every stochastic
replicate has more realized events. Every seed is retained, including
degenerate or failed replicates; do not replace seeds after seeing outcomes.

The divergent and baseline conditions share the same species/gene history
per seed. Likewise, divergent_turnover and turnover share history. Verify
matching T/G biological outputs rather than assuming equal seeds suffice;
parameter filenames or provenance differences are not biological differences.
This is 40 native simulation configurations and 30 derived datasets, for
70 condition/seed evaluations. Reuse immutable matching history checkpoints
only if their provenance and equivalence are verified.

### Missing Data

For missing20, rank all baseline extant gene IDs by SHA-256 of
`master_seed:missing20:gene_id`, breaking the negligible hash-tie case by ID.
Remove the first `floor(0.2 * N)` genes. The operation uses no gene-family
labels, sequences, truth pairs, or method outputs. Do not rescue deleted
orthologs or enforce favorable reference coverage. Project the true pair
set onto retained genes; record removed IDs and resulting eligible true pairs.
This represents missing observations, not evolutionary loss or fragmentation.

### Uneven Sampling And Count Control

Using the baseline species tree only, enumerate internal clades with two
to four extant descendants, excluding the root. Choose the largest such
clade; break ties by its sorted descendant-name tuple. Retain its
lexicographically first species and remove its other descendants. If no
eligible clade exists, record the condition as structurally inapplicable;
do not substitute a favorable clade. This retains five to seven taxa and
creates a localized sampling gap independently of method outcomes.

For taxon_count_control, remove the same number of taxa by ranking all
baseline species by SHA-256 of `master_seed:taxon_count_control:species`.
This is a count-matched hash sample, not guaranteed balanced sampling.
Report overlap with the clade-based removal. Restrict sequence inputs and
true pairs consistently and record every removed species/gene. Retained
species-tree branch lengths and any unary-node collapsing must be audited
if the tree is used as a supplied-tree diagnostic; primary methods infer
their own species trees.

## Frozen Methods And Endpoints

Run OrthoHMM high sensitivity and satellite_v2 from the same frozen
`7f3a9e4` source/configuration as prospective YGOB validation, and full
OrthoFinder 3.1.5 using DIAMOND. Extract its sequence-only MCL checkpoint
as a diagnostic. Use identical per-dataset FASTAs, four search CPUs and
four threads per OrthoHMM worker, and OrthoFinder `-t 4 -a 4 -S diamond`.
Record seed, package/tool versions, exact commands and failures. The CPU
budget is matched within this panel but is not matched to historical large
dataset timings. Resource measurements from contended runs are descriptive.

Primary targets are cross-species ortholog pairs whose true common ancestor
is a speciation event. Use event-derived pairs for phylogenetic methods;
expand high-sensitivity groups and the OrthoFinder checkpoint into
cross-species pairs. Reject unknown IDs; deduplicate pair orientation
explicitly. Score all predicted pairs within the retained input universe,
including between unrelated origin families. Such false positives must not
be projected away. Report pair micro F1/P/R, eligible true pairs, input and
prediction coverage, family copy counts, failures, wall time and memory.
Do not use root-origin family cliques as a resolved orthology or root-HOG
reference. Root-HOG scoring remains a separate, not-yet-defined endpoint.

Primary comparisons are the two OrthoHMM modes versus full OrthoFinder in
each of seven conditions. The sequence checkpoint is diagnostic only.
For each condition report all ten seed-level outcomes and their mean F1/P/R;
the mean across seeds is not the pooled-pair statistic. Estimate uncertainty
by paired resampling of whole seeds (20,000 PCG64 draws, seed 20261031),
recomputing mean differences. Report nominal 95% percentile intervals and
Bonferroni-adjusted F1 intervals across 14 planned comparisons. Precision
and recall intervals are exploratory and must be labeled as such. Ten
replicates yield limited tail resolution; do not claim exact simultaneous
coverage or calculate pair-independent significance tests.

If a method fails a condition/seed, report failure and do not manufacture
zero predictions as its accuracy result. Show failure fractions and
complete-case paired estimates with the excluded seed IDs; explicitly
state that these are conditional estimates and may be biased. Do not
silently reduce the requested panel or call it complete while failed
replicates remain unexplained. Full end-to-end reruns require documented
infrastructure justification and retain prior evidence.

## Launch And Completion Gates

Before launch, generate and commit a complete executable manifest with
expanded native parameters, source/environment hashes, exact commands and
immutable output destinations. Freeze and test missingness/taxon transforms,
pair scoring, and seed-level aggregation before reading method outcomes.
Verify extinct and single-survivor cases; observed missing survivor sequences
remain hard failures. No condition may be relaxed because a method scores
poorly.

The panel does not cover indels, domain shuffling, lineage-specific
composition models, transfer, species extinction, or structural fragments.
State these limitations. It also does not replace separate species-tree
error/parameter-neighborhood experiments, representative-size scaling,
independent curated validation, or the biological application. Completion
requires retained outputs, verified truth, all planned method evaluations,
failure accounting, paired summaries, and generated figures, not merely
parameter files or successful simulator exits.
