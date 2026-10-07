# SwissTrees Model-Based Divergence: Prospective Descriptive Protocol

## Scope And Selection

This implements the existing publication goal4.3 divergence/error-stratum
requirement, not a new accuracy campaign. All 18 retained SwissTrees reference
families and all 563 canonical members are selected, irrespective of prediction
outcomes. Their scores are already development-exposed. No default, threshold,
benchmark, alignment or orthology prediction is changed. The new feature stage
does not read prediction counts, labels, subgroup effects or native outputs.

Reuse the existing aligned.faa files, without new extraction, MAFFT, trimming,
masking, site selection, taxon filtering or reexecution of the old admission:

- corrected_swiss_identity_prepared_22102.json:
  SHA256 7bec0f40c8227f230bbf06cb52622947e256f2edfc2c88b991393fa570603e44.
- corrected_swiss_identity_admission_22102.json:
  SHA256 cb04162af62fbd58fcfe8f02cbb78bc53bcf49ac20a487aae89911bc3a4e7b2d.

Check these two direct reports, the 18 selected alignment bindings and their
admitted canonical memberships. Do not rehash the entire historical 194-file
admission or reread the original 78 proteomes. Alignment names must be unique;
each family must retain its full membership and admitted dimensions.

## Model And Execution

Use the already-installed static x86-64 IQ-TREE 3.0.1 binary:
`/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/iqtree-3.0.1-Linux-intel/bin/iqtree3`.
It is 11333032 bytes, SHA256
40424ccdb1d79c304641f910cb6c172ebb670214e50958352018e3ff9906ab8f.
No installation or update is authorized for this bounded analysis.

Prospectively use exactly one inference per family, in lexicographic family
order, with fixed WAG+G4, protein type, seed20261007, one thread, 4G application
memory ceiling and identical-sequence retention:

```text
iqtree3 -s ALIGNMENT --seqtype AA -m WAG+G4 --seed 20261007 -T 1
        --mem 4G -keep-ident --prefix FRESH_FAMILY_DIRECTORY/inference
```

All other search parameters remain IQ-TREE3.0.1 defaults. No ModelFinder,
topology constraints, species tree, supplied truth tree, support bootstrap,
branch support test, date calibration, relaxed clock, replicate selection or
model comparison. WAG uses its fixed model frequencies, not +F or +FO; G4 uses
four discrete gamma rate categories with the default mean approximation.

Use a new 2-CPU/8-GiB/4-hour, no-requeue Slurm job, sequential families and
600-second per-family timeout. Sanitize Python/dynamic-loader overrides and
limit OMP/BLAS threads to one. Record the source commit, exact source/protocol/
batch/binary bytes, Python/Biopython versions, job identity, commands, return
codes, timestamps and every family output/log. A fresh directory is required;
never overwrite or resume an existing attempt. Preserve any interrupted or
failed attempt and partial files. Continue to the other selected families after
a family failure, but do not retry or change that family's model/settings.
Ordinary shared-host contention is accepted, not a launch gate. Defer only for
unsafe memory or invalid allocation. This is descriptive feature-construction
cost, never OrthoHMM inference timing or isolated tool-speed evidence.

## Features And Bins

Parse each successful inference.treefile with Bio.Phylo. Require exactly the
admitted family leaf set, finite nonnegative lengths on every non-root edge
and finite nonnegative root length if present. Missing or invalid non-root
lengths, duplicate/missing/extra tips or multiple trees are failures. Do not
clip lengths. The root's incoming edge is excluded from pair distances.

For every unordered pair of retained family tips, calculate the patristic
distance (sum of edges on the tip-to-tip path), in model-estimated expected
amino-acid substitutions per site. Do not divide again by protein/alignment
length, normalize by a family maximum or call these biological time. Preserve
all pair values in a TSV. The primary family feature is the ordinary median
of ALL unordered pair distances, including within-species/paralog pairs and
zero distances. Also report the mean, min, max, pair count and tree length as
descriptors, not competing bin definitions or alternative accuracy endpoints.

Only if every selected family succeeds, compute the ordinary median of the
18 primary family features. Freeze bins: lower_or_equal_median = feature <=
this cutoff; higher_than_median = feature > this cutoff. Ties stay in the lower
bin, never split by outcomes. Retain an overall all-family bin. Empty bins
remain explicit with NA statistics; no outcome-driven cutoff alternatives.
If any family fails, retain its missing feature and failure, do not compute
a successful-only cutoff or project subset scores. No automatic retry.

An independently implemented reader must verify selected input/source/output
bindings, full tip sets, every pair distance by edge-split summation (not the
exporter's path-distance implementation), all descriptors and the exact
median/tie bins. Shared Bio.Phylo parsing is disclosed; this is independent
arithmetic, not independent phylogenetic inference/model/alignment validation.

## Subsequent Accuracy Projection

Only after feature readback succeeds may a separate stage read the existing
native_qfo_three_cell_strata_20261007_v1/report.json (SHA256
55ed72ce771c756424877cfa412da329ff104efb7dca4b35113b0d80703a59e5)
and its v2 independent readback (SHA256
6f3eef63c0e96e1d85805fcfb70ed9b66fbb1db76a08ca2a105483422a0f6969).
Use all 54 retained family count rows, the exact same 18 family/member sets,
the already-defined TP/2+1, FP/2+1, FN/2+1 prior and macro-family precision/
recall followed by harmonic F1. Preserve TN; do not pool pairs or average
family F1. Export three cells x three bins and both conditional contrasts
(R_at_P0_C0 and C_at_P0_R0), with differences in percentage points. No new
bootstrap draws, intervals, significance claims, C-by-R interaction claim,
timing admission or reexecution of a raw scorer. Independent rational
readback must cover every score/contrast row. Cell7's failed timing remains
ineligible. Keep negative and neutral effects and all family descriptors.

## Limitations And Primary Documentation

These are estimated, alignment-dependent, fixed-model-dependent distances on
exposed families, not known history, validated elapsed evolutionary time,
causal explanations, model-fit confirmation or generalization evidence.
Within-family divergence combines speciation, paralogy, taxon sampling,
domain/length differences, compositional effects, gaps and alignment error.
Fixed WAG stationarity/homogeneity assumptions and heuristic tree-search
uncertainty remain; no support intervals or model-adequacy tests are inferred
from successful execution. Alignment-quality warnings are retained, not used
to exclude families after seeing results. All current descriptor and original
TreeFam/reference limitations remain explicit.

Primary documentation inspected2026-10-07:

- [FastTree author documentation](https://morgannprice.github.io/fasttree/):
  CAT-derived Gamma20 estimates have small-alignment limitations. Fourteen of
  these 18 retained families have fewer than 50 sequences. This motivates
  using the existing IQ-TREE installation rather than FastTree CAT rescaling.
- [IQ-TREE command reference](https://iqtree.github.io/doc/Command-Reference):
  explicit model, seed, threads, memory, prefix and identical-tip retention;
  interrupted commands can auto-resume, hence the strict fresh-directory rule.
- [IQ-TREE substitution models](https://iqtree.github.io/doc/Substitution-Models):
  empirical protein matrices/frequencies and discrete gamma rate categories.

Commit/push this protocol before selected inference or feature computation;
commit/push tested new sources before executing them. No existing frozen
scientific source, native23902 or terminal-review23910 is changed by this work.
