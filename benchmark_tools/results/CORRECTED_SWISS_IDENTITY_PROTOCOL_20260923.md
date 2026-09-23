# Corrected SwissTrees Identity Strata

Exploratory, development-exposed analysis specified after aggregate outcomes
have been inspected. No method tuning or new independent-validation claim.
Freeze this plan and implementation before generating alignments or joining
these features to prediction errors. Preserve all 18 families and failures.

Use the corrected sequence inventory SHA-256
912363276fa5456f11aeff9f398fce71aec8d377c8c8eda7b144105945e554f1,
its exact 78 FASTA identities, and its 18 reference-family memberships
(563 distinct proteins). Read no method predictions or confusion counts.
Require exact accession coverage and unchanged descriptor lengths. Reuse
the validated OrthoBench MAFFT helper, SHA-256
9a787bffec2a9481c8987a7e7c5066982c1ec3fe4c3f571cf7a02e69162e885d,
and the frozen MAFFT 7.525 installation. Record its complete helper-binary
inventory before and after execution.

For each family sort accessions; uppercase proteins, remove explicit stop
symbols with counts retained, preserve X and reject other noncanonical
symbols. Do not modify inference FASTAs. Run MAFFT --amino --auto --thread 1,
with at most four families concurrently. Require alignment membership and
ungapped normalized residues to match the exact inputs. Preserve failures
without automatic retry; no partial successful subset becomes an admitted
feature panel. Scheduler allocation is 4 CPUs, 16 GiB, two hours, no requeue.
Preparation timing is descriptive, not comparative inference timing.

For every unordered protein pair, identity is identical canonical residues
divided by aligned positions where both residues are canonical. Exclude gaps
and X from both numerator and denominator. Any pair with zero comparable
positions makes the whole family's identity missing. Otherwise use the
unweighted arithmetic mean of pair identities. Freeze the median of family
means as the split; ties are lower-identity. Retain lower, higher and missing
bins, including empty bins. Sequence identity is alignment-dependent and
not a calibrated evolutionary distance or evidence of duplication history.

After complete input-only preparation and admission, display all admitted
corrected methods in each bin, including missing OrthoMCL results. Use the
existing raw/2+1 per-family SwissTrees precision/recall convention, arithmetic
macro P/R across selected families, and their harmonic F1. Do not pool pairs
or average family F1. Report family membership, number of families and
descriptive differences against full OrthoFinder. No additional bootstrap,
p-value, significant-subgroup claim or changed multiplicity denominator is
planned for these exploratory displays. Existing primary composition and
overall comparator uncertainty analyses remain separate and unchanged.

These bins can differ in taxon coverage, size, domains and composition.
Associations do not identify mechanisms. This does not supply independent
fragment annotations, true divergence estimates, or complete the broader
error-analysis requirement by itself. No outcomes may be used to revise bins.
