# Recovered SwissTrees Uncertainty Protocol

Freeze this analysis before calculating or inspecting bootstrap intervals.
The four stage point estimates have already been inspected and reported;
this is development-exposed follow-up analysis, not prospective independent
confirmation. No stage or family is chosen from its observed effect.

## Verified Sufficient Statistics

Input qfo_swiss_counts_20260917.json SHA256:
546bb5bd6957c8ea990324b79fa31f22b0ed721bc7d6b94b609ab15258f97183.
The original admission SHA256 is
89b683d0fc7fe9964ce5b6182832bb8eb28758fad3f6bbb4d706926143805956.

All18 reference families and10765 relation labels occur in every stage.
Reference truth and represented gene inventories are identical across stages;
there are no shared represented genes across families. The frozen Darwin
runtime independently loads the frozen reference and confirms that every
relation is stored in increasing protein-ID orientation only. Each family has
more than5 mapped proteins and raw records cover every mapped protein and
relation. Duplicate/reversed duplicate raw rows are rejected.

The actual container RefPhyloTest.drw matches the inspected repository scorer
SHA256262c8d1f06527461e11391682ccaccc637675f84a5064c1658b80bf1b956bd51.
It adds length(matches)/2 and starts each confusion count at1. Therefore,
for this one-direction reference, the native count is raw_count/2+1,
equivalently raw_count+2 for computing ratios. A first audit using raw_count+1
failed to reproduce the native metrics; that interpretation was rejected,
not used to alter the benchmark. The corrected reconstruction reproduces
all family and aggregate native coordinates within5e-8 absolute tolerance.

For each family compute precision TP/(TP+FP), recall TP/(TP+FN) using the
native pseudocounts. The benchmark aggregate is the unweighted mean precision
and unweighted mean recall across families. Its project-reported F1 is the
harmonic mean of those two means, NOT mean family F1 or a pooled-pair F1.

## Fixed Resampling

- Units: all18 reference families, sampled with replacement; same family
  multiplicities for every stage and contrast. Never resample individual pairs.
-100000 draws, numpy.Generator(numpy.PCG64(20260919)), multinomial counts
  with18 trials and equal probability1/18 for each family.
- Recompute mean precision and recall under every draw and then their harmonic
  mean. Use full-precision count-derived statistics, retaining the small native
  output-rounding differences rather than copying rounded aggregate scores.
- Fixed contrasts, candidate minus reference: strict_profiles minus multipass;
  strict_profiles_refined minus multipass_refined; multipass_refined minus
  multipass; strict_profiles_refined minus strict_profiles.
- Report F1, precision and recall for each contrast in raw0-to1 units.
  Report nominal95% paired percentile intervals and Bonferroni-adjusted
  percentile intervals for the12 planned contrast/metric endpoints, with
  quantiles0.05/(2*12) and1-0.05/(2*12). Use numpy linear quantile interpolation.
- Report all family effect records and descriptive wins/ties/losses, defining
  ties as absolute difference<=1e-10. Do not add significance claims for these
  descriptive counts or drop neutral/negative cases.
- Input count report must pass its pinned hash. Reject missing stages/families,
  non-integer or negative raw counts, changed family membership, shared genes,
  differing reference truth totals, or disagreement with stored statistics.
  Failed validation produces no interval result or score substitution.

## Interpretation Limits

The12-endpoint adjustment covers only this fixed SwissTrees stage analysis,
not all publication experiments or previous parameter searches. With18 curated
families, exchangeability is approximate: shared evolutionary history and
merged predicted groups can correlate disjoint families. Percentile intervals
are approximate and not independent confirmation or model-selection-adjusted
evidence. Native reported stderr fields are not used as paired intervals.

Profiles include downstream graph and singleton-assignment changes; refined
means sequence-based post-clustering refinement, not phylogenetic inference.
No method retuning or claim of total HMM contribution is authorized. Other QfO
challenges and the secondary six-metric mean need separate defensible uncertainty
methods; SwissTrees intervals must not be propagated to them.
