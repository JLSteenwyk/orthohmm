# Native SwissTrees Duplication-Annotation Projection

Retrospective descriptive extension of goal4.3 to the two completed native
P0C0 cells. Initial HMM search remains on; profiles/candidate expansion off.
R0 predictions are group-clique pairs; R1 predictions are resolved native
pairs. This is not total HMM contribution, method tuning or a new independent
test. Commit this protocol before calculating the native bin outcomes.

Reuse the September23 feature definition and dual-traversal admission:
`swiss_duplication_features_v2_20260923.json`, SHA256
`97b0c4755d6a9df258d5c3f60fc0d5d25f1e5c09c42216c754a245a67d1942ec`.
Recompute count-derived exact rational fractions, median and tie-preserving
memberships. Retain all/lower/upper/missing bins, including empty missing.
The median stays7/48; lower includes APP/ASTER/BAMBI/BAR/CITE/Clusterin/NOX/
POP/SUMF; upper CASP/GH14/HOX/MAPT/PSEN/RPS/SERC/TRFE/VATB. No cutoff search.

Join only the completed native SwissTrees count bindings in the independently
read-back native sequence projection (SHA256
`4650f614d2fdeebd65cbcbf4999c61c83cc549560d2738bfbb5502f686f660e3`;
reader `8d0356bc0e159d204939271578e23a6a5b6efcf6198b75a18c0f3d38e0e67c0a`).
Require exact18-family/563-protein agreement with the retained mapping
(`0d3e736a350782609c68764bc19387500ea64bd5045c5fcfb41930dd89c1ce9d`).
Re-read both original native raw files and require identical reference truth,
complete memberships and all36 original family counts/statistics. Recompute
the existing prior, macro precision/recall and their harmonic F1, not pooled
pair F1 or arithmetic mean family F1. Emit eight rows/four R1-minus-R0
differences; empty statistics/differences null and TSV NA, not zero.

Primary Fraction/median implementation and separate stdlib reader using integer
cross-products for bin decisions and alternative raw-count/statistic arithmetic
must agree. Check source/binding/raw/output identities before and after reads.
Test ties, missing/empty bins, alias-sensitive informative-node denominator,
malformed counts/fractions, membership changes, row completeness, semantics,
changed scores/differences and output refusal. Preserve failed attempts.

Do not repeat original Darwin extraction or whole-study admission. Explicitly
inherit the old tree/mapping/runtime acquisition and traversal checks, rather
than claiming them freshly validated. The current mapping report and extractor
source are directly checked, not the full nested dependency graph. No historical
intervals transferred, new bootstrap draws, pair-IID inference or causal claim.
This reference-derived annotation fraction is related to the reference pair
labels; default-S is not explicit speciation, informative nodes need not equal
unique genes minus one, and annotation density is not a true duplication rate.
Development exposure, size/composition confounding and biological limits remain.

Use retained runtime without installs. Record shared-host postprocessing costs,
not inference costs or isolated speed. CPU/memory-bandwidth/I/O contention has
unknown, potentially tool-dependent effects. Failed R1 timing stays ineligible.
Leave original jobs22444/22445/22450/22451/22452 and frozen protocol bytes alone.
No publication-complete status from this bounded extension.
