# SwissTrees Domain-Strata Protocol

Freeze before calculating or inspecting domain-stratified accuracy contrasts.
Overall and family accuracy have previously been inspected in this project;
this is a retrospective, development-exposed explanatory analysis, not a
prospective independent test. No method retuning or family selection follows.

## Annotation Evidence

`swiss_domain_annotation_inventory_20260917.json` SHA256:
d5e269158c2fb603acc0805a1c140a88758342d7b8e32b3094750ade22fbc06c.
`qfo_swiss_comparator_counts_20260917.json` SHA256:
2995868b0407ceda2e7422db6c3fc3b99523716e4bec8b1b43e15f296e40769b.

The input-only inventory scans all 78 retained QfO FAS annotation JSON files,
pinning each file's bytes/hash, and extracts exact accession matches for all
563 represented proteins in the 18 SwissTrees families. All have annotations;
five have no Pfam hits, 225 have multiple distinct Pfam types and 63 have
multiple instances of at least one Pfam type. No fuzzy ID mapping was used.
Eight extraction tests cover repeat/type distinctions, missing versus zero-hit
annotations, missing namespace, invalid coordinates and duplicate instances.

This extraction uses reference membership but no prediction counts or scores.
Pfam hits provide direct annotation features rather than length-only proxies.
They are not validated domain-loss, fragment or ancestral duplication labels.
Coordinates are retained without calculating widths or sequence coverage because
coordinate conventions have not been independently verified. The same resource
underlies FAS, so do not claim independent validation of the FAS endpoint.

## Fixed Strata

Primary split, evaluated for all eight retained comparison methods:

- `median_pfam_types_below_two`: family median distinct Pfam types < 2.
- `median_pfam_types_at_least_two`: family median distinct Pfam types >= 2.

The annotation-only inventory gives 12 and 6 families respectively. The higher
stratum contains APP, BAR, HOX, NOX, TRFE and VATB. Retain every family, including
zero-hit protein records. Both bins have complete annotation coverage; abort
if the frozen coverage or family inventory changes rather than silently imputing.

Secondary descriptive split:

- `repeated_type_fraction_below_quarter`: fewer than 25% of annotated proteins
  have multiple instances of a Pfam type.
- `repeated_type_fraction_at_least_quarter`: at least 25% do.

The second bin contains MAPT, PSEN and TRFE; the other contains 15 families.
Report both bins, all family records and point estimates, but no interval or
significance claim for this three-family repeat subset. Repeated annotation
instances need not represent biologically validated tandem repeats.

## Statistics And Uncertainty

For each bin, recompute native per-family P/R from raw/2+1 confusion counts,
then unweighted macro P/R and their harmonic F1. Never pool pairs or average
family F1 as a substitute for this statistic. Display all eight methods.

Fixed inferential contrasts for the primary split, candidate minus reference:

1. High-sensitivity OrthoHMM minus full OrthoFinder.
2. Phylogenetic OrthoHMM minus full OrthoFinder.
3. Phylogenetic OrthoHMM minus high-sensitivity OrthoHMM.

Use 100,000 draws, NumPy PCG64 seed20260921. In the stated primary-bin order,
sample n families with replacement from each bin using multinomial(n,1/n)
multiplicities. Use the same draws for all methods and contrasts within a bin.
Recompute the actual macro/harmonic statistic in every draw. Also compute the
interaction as (candidate-reference in higher-type bin) minus
(candidate-reference in lower-type bin), pairing the corresponding independent
bin draws. This is a difference in contrasts, not a causal effect of domains.

Report F1, PPV and TPR differences with nominal95% percentile intervals and
Bonferroni intervals across27endpoints:3contrasts * 3metrics * (2bins+1interaction).
Use linear interpolation and adjusted quantiles0.05/54 and1-0.05/54. Retain all
intervals, including negative and zero-crossing results. Descriptive family
wins/ties/losses use absolute difference<=1e-10 as a tie; do not assign separate
inferential significance to those counts. Do not enlarge or replace strata
after inspecting outcomes.

## Interpretation Limits

Six and twelve curated families support only approximate conditional bootstrap
sensitivity estimates. Shared history and merged predictions may correlate
families. The annotation bins can also differ in divergence, taxon sampling,
family size and curation; stratification cannot establish the cause of an error.
The three-family repeat bin is descriptive only. Inclusion of zero is not
equivalence. Adjustment covers these27endpoints, not prior development or all
publication analyses. Preserve the MCL-checkpoint and supplied-tree FastOMA
diagnostic labels. The OrthoHMM configuration contrast changes more than
reconciliation. These estimates do not cover other QfO metrics or the secondary
six-metric mean, and do not fill the independent fragment/duplication-label gaps.
