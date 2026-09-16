# OrthoBench Error Analysis: Exploratory Protocol

Specified after development, independent YGOB evaluation and the completed
OrthoBench factorial/controls. This is not preregistered confirmation. Freeze
this feature/contrast plan before joining the new features to method errors;
do not modify the inference method or thresholds from these outcomes.

## Input Feature Preparation

Use all251,378 proteins from the twelve hash-pinned OrthoBench proteomes.
The first-stage preparer reads only the input-manifest FASTA records, not
predictions, benchmark scores or reference-family membership. Retain one TSV
row per unique protein and exact input/output hashes. Refuse duplicates,
changed inputs and an incomplete universe; preserve unsuccessful output
without a completion manifest.

Residue length counts ASCII letters A-Z after case normalization. Separately
record canonical20-residue count, noncanonical letters, stops, gaps and other
symbols. Normalize Shannon entropy of canonical-residue frequencies by
log2(20); no canonical residues means missing, not zero entropy. Frequencies
are calculated in a fixed residue order for reproducibility.

Two descriptive flags are frozen: short sequence means residue length<100
(missing for zero residues); composition-concentrated means normalized
entropy<0.8, evaluated only with at least20 canonical residues and canonical
fraction>=0.9. Otherwise composition status is missing. These are exploratory
descriptors, not validated fragment/low-complexity/domain classifiers.

## Planned Family Strata

Join the verified input features to the frozen70 RefOG membership records,
without using any method's predictions to define strata. Report every family
and its features before producing stratified outcome tables.

- Family size:2-20,21-50,>50 members, retaining the existing diagnostic's bins.
- Copy number: at most one versus more than one reference member per species.
  This is a copy-number descriptor, not an inferred duplication history.
- Alignment identity: lower/higher mean pairwise identity, split at the
  median of available family values (ties in lower), plus a missing category.
  Validate or rebuild reference-family alignments against the input sequences
  using pinned MAFFT. Calculate each pair's identity only at positions where
  both residues are canonical; a pair with no comparable positions makes the
  family statistic missing. Record alignments and their provenance. Do not
  reuse older identity summaries without these checks.
- Relative length: any member below0.5 times its family's median residue
  length versus no such member, plus missing for any zero-residue member.
  This is a short-relative-to-family descriptor, not evidence of fragmentation.
- Composition: any member with the frozen concentrated-composition flag versus
  all members with an evaluated false flag, plus missing when none is true
  and one or more members are unevaluable.

These five dimensions have3+2+3+3+3=14 planned strata. They overlap; neither
their sample sizes nor their effects may be pooled as independent evidence.
Actual duplication history, annotated fragments, domain architecture and
biological mechanisms still need separate supporting data. These proxies do
not satisfy those requirements by themselves. QfO extension also remains open.

## Outcome Comparisons

Use the already frozen publication outputs for OrthoHMM high-sensitivity,
OrthoHMM satellite_v2 phylogenetics and full OrthoFinder3.1.5. Native semantics,
input coverage and conversion checks must pass before outcome assembly.
Within each stratum, report both OrthoHMM configurations minus full OrthoFinder
for F1, precision and recall. Keep the original full-reference per-family
TP/FP/FN and low-certainty conventions, then recompute the weighted statistic
from the selected family records. Do not silently alter the reference universe
or count only pairs between selected-family genes.

Use20,000 paired RefOG bootstrap replicates, seed20260918, within each stratum.
Report nominal intervals and Bonferroni tail adjustment over all84 planned
endpoints (14 strata x2 contrasts x3 metrics), retaining that denominator
when a stratum is empty, missing or nonestimable. Bins with fewer than five
families receive descriptive estimates only, no bootstrap claims. This rule
does not imply adequate power in larger bins. Empty bins remain explicit.
Report family counts, coverage, weighted counts and family-level wins/ties/
losses; avoid gene-pair independence assumptions and equivalence claims.

## Mechanistic Tracing

Trace all70 reference families through retained initial search, candidate
membership and final grouping/reconciliation where records permit. Record
missing intermediate evidence rather than assuming a mechanism from score
changes. For concise illustrations, use the smallest SHA256-ranked RefOG ID
within each occupied feature stratum, deduplicating shared selections. Keep
the full trace table available, including neutral and adverse cases. These
development-exposed traces are not a separate biological application and
do not use competitor predictions as ground truth.

## Sequence Feature Preparation Completed

Frozen executor8bccb33, job21308, completed0:0 in23 seconds. All251,378 unique
proteins are represented. Source/input/table hashes and per-proteome row
counts were independently verified. The32,500,737-byte TSV is retained outside
Git with SHA256d315da458240fe2020b128b606adfcfdf15d4de34ee937a6b2aa7735f7452f97;
`ob_sequence_features_prepared_20260916.json` records provenance and counts.

Descriptive input counts:13,822 sequences below100 residues;2,345 with the
composition-concentration flag and460 unevaluable for that flag. There are
3,170 proteins with noncanonical letters and177 with stop symbols; no empty
residue sequences, gaps or other symbols were observed. No family-stratum
outcomes, biological error explanations or accuracy effects are inferred
from these counts. Reference-alignment admission and family-feature joins
remain to be completed before the prespecified outcome analysis.
