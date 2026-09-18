# Corrected SwissTrees Sequence Strata

Freeze before calculating corrected-release composition-stratified outcomes.
Historical and some corrected aggregate scores have already been inspected;
this is an exploratory development-exposed analysis, not independent testing.
Do not tune the method or choose bins from these outcomes.

## Input-Only Preparation

Read the78 corrected staged FASTAs with manifest SHA-256
`07a890eb816944f946d46039559a6046c2a3664f033eaa0f30f35bb25b9c9ab8`.
Use only the18-family,563-protein membership from the existing SwissTrees
count inventory2995868b0407ceda2e7422db6c3fc3b99523716e4bec8b1b43e15f296e40769b.
Do not use historical counts as corrected outcomes. Require exact accession
matches and complete coverage, check FASTA hashes before/after extraction,
and preserve the549 previously matched descriptors byte-for-value.

Reuse the frozen canonical-residue entropy and literal fragment-label
extractor. Global Shannon entropy excludes noncanonical characters and is
normalized by log2(20). Length includes all sequence characters as in that
extractor. Retain descriptions, noncanonical counts and lengths. A protein's
entropy is eligible for stratification only with at least20 canonical
residues and canonical fraction>=0.9. A family with any ineligible protein
has missing family entropy; never silently drop members.

Primary split: compute each fully eligible family's median normalized
entropy, then split at the median of available family medians. Ties go to
the lower bin. Preserve an explicit missing bin. No outcome information
determines this cutoff. Freeze actual memberships in the resulting inventory.

Secondary descriptive bins: any eligible member entropy<0.8 versus all
evaluated false, retaining missing when none is true and some are ineligible;
any member length<0.5 times family median versus none; any literal parenthesized
Fragment/Fragments description versus none. These overlapping bins receive
no additional inferential claims. Absent fragment text does not prove complete
sequence, and relative length does not diagnose fragmentation. See the
[UniProt manual](https://web.expasy.org/docs/userman.html) for annotation flags;
the retained FASTA descriptions, not current UniProt entries, are the data.

## Planned Outcomes

Wait for admitted corrected-release predictions/counts. Display all available
methods with honest missingness, but do not substitute original-release scores.
Primary contrasts: high-sensitivity OrthoHMM minus full OrthoFinder,
phylogenetic OrthoHMM minus full OrthoFinder, and phylogenetic minus
high-sensitivity OrthoHMM. The configuration contrast is not a pure ablation.

Within each primary bin, recompute native SwissTrees per-family smoothed P/R
from the frozen scoring convention, then macro P/R and harmonic F1. Do not
pool pairs or average per-family F1 instead. With at least five families in
each bin, use100000 PCG64 multinomial paired family draws, seed20260924,
lower bin then higher bin, same draws across methods within each bin.
Report lower/higher-bin contrasts and higher-minus-lower interactions for
F1,PPV,TPR. Use nominal95% and Bonferroni27 intervals, linear quantiles
0.05/54 and1-0.05/54. Retain the27-endpoint adjustment if any endpoint is
nonestimable. Any bin with fewer than five families is descriptive only;
an interaction requires both bins eligible. Missing-entropy and all secondary
bins are descriptive only. Retain all directions and zero-crossing results.

These small curated-family bootstraps are conditional and approximate;
shared history/predictions can correlate families. Composition strata may
also differ in taxa, domain architecture, size and divergence. Associations
do not establish mechanisms or performance on new families. Adjustment does
not cover prior development or all publication analyses. No other QfO
endpoint or the secondary mean inherits these uncertainty estimates.
