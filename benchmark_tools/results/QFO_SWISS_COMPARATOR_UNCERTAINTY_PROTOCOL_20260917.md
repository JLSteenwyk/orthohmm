# SwissTrees Comparator Uncertainty Protocol

Freeze before calculating or inspecting intervals. Historical point estimates
have already been inspected and used during development. This is retrospective,
development-exposed inference, not independent confirmation or adjustment for
past model selection. Retain every comparison regardless of direction.

## Verified Inputs

`qfo_swiss_comparator_counts_20260917.json` SHA256:
2995868b0407ceda2e7422db6c3fc3b99523716e4bec8b1b43e15f296e40769b.
The parent eight-method comparison SHA256 is
094f842ad8211450519d4238c0b3a5020a7969049d4b7306f5465d5acbf914a3.
The independent native/reference count audit SHA256 is
546bb5bd6957c8ea990324b79fa31f22b0ed721bc7d6b94b609ab15258f97183.

All eight retained methods have the same 18 reference families, represented
members and 10,765 pair truth labels. Represented genes do not overlap between
families. Raw counts reproduce every retained aggregate precision, recall and
harmonic F1 within 5e-8. The reference audit established increasing-ID relation
orientation and native per-family confusion counts of raw_count/2+1.

The aggregate is mean family precision and mean family recall, followed by their
harmonic mean. It is neither pooled-pair F1 nor the average of family F1 values.
Use the full-precision reconstructed statistic in every replicate.

## Fixed Contrasts

Candidate minus reference, using these eight contrasts:

1. orthohmm_high_sensitivity minus orthofinder_3_1_5_full.
2. orthohmm_phylogeny_satellite_v2 minus orthofinder_3_1_5_full.
3. orthofinder_3_1_5_sequence_only minus orthofinder_3_1_5_full.
4. sonicparanoid_2_0_9 minus orthofinder_3_1_5_full.
5. proteinortho_6_3_6 minus orthofinder_3_1_5_full.
6. fastoma_0_3_5 minus orthofinder_3_1_5_full.
7. orthomcl_1_4 minus orthofinder_3_1_5_full.
8. orthohmm_phylogeny_satellite_v2 minus orthohmm_high_sensitivity.

The third is an MCL-checkpoint diagnostic, not a separate final OrthoFinder
output. FastOMA used a supplied OrthoFinder tree. The eighth is a configuration
comparison, not a pure reconciliation ablation. Preserve these distinctions.

## Fixed Resampling

- Sample all 18 families with replacement using shared multiplicities across
  every method and contrast. Never treat dependent gene pairs as independent.
- Use 100,000 draws from numpy.Generator(numpy.PCG64(20260920)), each a multinomial
  sample with 18 trials and probability 1/18 per family.
- Recompute weighted mean precision/recall and their harmonic mean in every draw.
- For F1, PPV and TPR report observed differences in raw 0-to-1 units, nominal
  95% paired percentile intervals and Bonferroni-adjusted percentile intervals
  across all 8*3=24 endpoints. Adjusted quantiles are 0.05/48 and 1-0.05/48;
  use NumPy linear interpolation.
- Report every family's F1/PPV/TPR differences and descriptive wins/ties/losses;
  ties are absolute differences at most 1e-10. Do not select families or assign
  inferential significance to descriptive counts.
- Validate pinned input identities, exact method/family inventories, finite
  nonnegative integer counts, shared truth totals and member sets, disjoint
  represented families, and stored family/aggregate statistics. Fail closed
  on missing or inconsistent evidence; do not substitute rounded scores.

## Interpretation

Only 18 curated families are available. Disjoint membership does not establish
biological independence: shared evolutionary history and merged predictions can
correlate families. This family bootstrap is an approximate conditional
sensitivity analysis. It does not prove generalization to unrelated families
or account for prior tuning. State these limitations beside the results.

Multiplicity adjustment covers these 24 endpoints, not every publication
analysis. Do not combine with the separate four-stage interval analysis as
independent evidence. These historical prediction outputs differ from recovered
ablations; do not transfer scores or intervals between them. No conclusions
about other QfO challenges, the secondary six-metric mean, global superiority,
or equivalence follow from these intervals. No method retuning is authorized.
