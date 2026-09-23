# Retained-Tree Duplication Feature Protocol

Freeze this definition before extracting features or joining prediction
outcomes. This is a descriptive extension on development-exposed QfO, not
a new confirmatory hypothesis or independent test set.

## Input Gate

Require exact reconstruction of all 18 families from the retained 2020
reference and identifier mapping. Pin `swiss_retained_mapping_v2_20260923.json`
to SHA256 `0d3e736a350782609c68764bc19387500ea64bd5045c5fcfb41930dd89c1ce9d`
and its reconstruction helper to
`871acee01b9e21441483230076cdadbc4d0d76df771f0fcfde5e3168683029ae`.
Recheck the input identities and recompute native traversal/mapping; require
the family diagnostics to agree before admitting features. Do not substitute
current SwissTree downloads, inferred tool trees or selected favorable families.

## Feature

Traverse the retained binary tree in its stored orientation. Map leaves by
the validated exact-identifier/first-underscore-alias rule. At each internal
node, obtain the mapped descendant sets of each child. Preserve the native
left-minus-right rule for overlapping sets. Count a node as informative only
when both resulting sets are nonempty. This excludes events supported solely
by taxa absent from the benchmark and excludes zero-cross-pair nodes.

Classify informative nodes using the native annotation precedence:
explicit speciation substring before duplication substring within each
annotation, annotations in order; numeric-only annotations are skipped.
No event annotation defaults to S for reference reconstruction, but record
this separately as `default_speciation_nodes`. Record explicit S and D counts,
informative-node count and child-overlap count. Unsupported annotations fail
closed. Preserve compound and wrapped NHX forms.

The sole primary feature is
`duplication_fraction = explicit_duplication_nodes / informative_nodes`.
Use exact rational arithmetic for bin assignment. The denominator is the
observed informative-node count, not `n_genes - 1`, because identifier
collisions can make an ordinary unique-leaf-tree assumption inappropriate.
This is an annotation fraction, not a duplication rate per evolutionary time,
a count of all ancestral duplications, or a mechanistic explanation of errors.

## Strata and Missingness

Take the median of the available family fractions without reading prediction
results. Lower stratum: fraction at or below the median. Upper stratum:
strictly above it. Never split tied values to force equal bin sizes. Retain a
missing bin for zero informative nodes and report the full 18-family result.
An invalid mapping or failed reconstruction blocks feature admission rather
than turning that family into missing. Report memberships, exact median,
raw counts and fractions. No additional cutoff search or exclusion analysis.

## Scores and Interpretation

Only after feature admission, join the frozen corrected seven-method counts,
retaining unadmitted OrthoMCL as NA. Recompute native per-family priors,
macro precision and macro recall, then their harmonic mean for F1. Use full
OrthoFinder as the difference reference. Retain all eight method rows in all,
lower, upper and missing bins, even when a bin is empty. Record each method's
prediction semantics. Do not use arithmetic mean family F1 or pooled pairs.

Report descriptive F1, precision and recall and their differences; no new
p-values, bootstrap intervals or subgroup superiority claims. Explain the
relationship between this reference-derived feature and reference pair
labels, possible family-size/composition confounding, alias collisions and
default-S annotations. Do not claim this establishes causation or generalizes
outside these families. Existing global uncertainty analyses remain unchanged.
