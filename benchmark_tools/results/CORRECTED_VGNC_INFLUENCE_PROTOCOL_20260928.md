# Corrected VGNC Deletion Sensitivity

Exploratory extension of the historical four-stage deletion diagnostic, defined
before this corrected-panel calculation. Retain every one of the 16,844
reference overlap blocks and all eight corrected methods. For each method,
delete each block separately, removing every incident scored cell once for
that deletion (including cross-block FP cells). Recompute precision, recall
and F1 from remaining counts. Do not rerun native inference or eligibility.

Compare every nonbaseline method with full OrthoFinder on the same deleted
block. Retain all deletion rows, extrema and positive/zero/negative difference
counts; report all seven contrasts, not a subset selected by outcome. Preserve
full original counts and scores. Validate incident counting on synthetic cases
and independently recompute native deletion contrasts from sparse cells.

These ranges are not confidence intervals. Deletions are dependent and change
the scored table. Single-block stability does not imply robustness to joint
deletions, uncertainty coverage, new families or new clades. No defaults,
benchmark endpoints or native VGNC uncertainty estimates are changed.
