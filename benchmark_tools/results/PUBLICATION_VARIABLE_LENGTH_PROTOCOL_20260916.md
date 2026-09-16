# Prospective Heterogeneous-Family-Length Simulation Panel

## Motivation and Scope

The original fixed-300-residue panel remains frozen and is retained as a
stress test. Native OrthoFinder 3.1.5 normalization overflow was reproduced
before accuracy scoring; see `SIMULATION_NATIVE_FAILURE_AUDIT_20260916.md`.
This additional panel is motivated by that failure diagnosis, not a search
for parameters that favor either method. Do not pool its results with the
original panel or describe the original panel as successful validation.

The additional panel introduces between-family length heterogeneity while
retaining the original evolutionary conditions, WAG model, native tree and
gene-history simulators, four method configurations, CPU allocations, truth
definition and conversion semantics. It is a simplified simulation, not an
empirically fitted distribution of biological protein lengths. It still
does not model indels, within-family length changes, domains, composition
bias, transfer, or species extinction.

## Frozen Design

- Ten new seeds: 20261101 through 20261110 inclusive.
- Eight extant species and 100 root-origin families per native simulation.
- The four native conditions retain the original D/L and sequence scaling:
  baseline (2,1,0.2), divergent (2,1,0.8), turnover (10,8,0.2), and
  divergent_turnover (10,8,0.8). Other T/G/S settings are unchanged.
- For each seed and canonical positive-integer family ID, encode the ASCII
  string `orthohmm-sim-length-v2:{seed}:{family}`; take SHA-256, interpret its
  first eight bytes as an unsigned big-endian integer, and set amino-acid
  length to `100 + integer % 401`. Thus lengths lie in [100,500], with an
  approximately uniform hash allocation centered on 300. Do not redraw or
  reject lengths based on tree histories, tool completion, or scores.
- Assign one length per root-origin family, shared by all its descendants,
  paralogs, and the four native conditions for that seed. No family
  origination is enabled. Materialize and hash each mapping before generation.
- Retain the three baseline-derived conditions: missing20, uneven_taxa and
  taxon_count_control, with the original deterministic rules and the new seed.
  This gives 40 native runs and 70 datasets, without replacing failed seeds.
- Require identical T/G biological histories for the same 20 within-seed
  baseline/divergent and turnover/divergent_turnover comparisons, plus equal
  family-length mappings. Verify every exported sequence against its assigned
  family length and retain the existing event/XML/tree/FASTA truth checks.

## Implementation and Validation

Pinned Zombi's native Sf mode uses a codon model, including when amino-acid
sequence mode was requested. Therefore use ordinary S mode with a narrowly
scoped wrapper: set `SequenceSimulator.size` from the frozen family mapping
for each call to its original `run` method, then restore the previous value.
Do not change Pyvolve scoring/evolution, WAG parameters, tree branch lengths,
history files, or native name correction. Keep the existing explicit family
random-seed adapter. Default fixed-length execution remains unchanged.

Before generating this panel, use seed 20260918 for an eight-species,
100-family baseline smoke test. It is not a scientific-panel member. Require
repeat-run identity, heterogeneous exported lengths matching the mapping,
truth consistency, and native completion/numerical checks for the existing
OrthoFinder command. A failure is diagnostic evidence, not permission to
choose a different length rule without a documented new protocol. Do not
calculate smoke accuracy or select the panel based on an accuracy outcome.

## Analysis and Claim Boundaries

Keep OrthoHMM core revision 7f3a9e4 and the original high-sensitivity and
satellite_v2 settings. Keep unmodified OrthoFinder 3.1.5 full inference and
its sequence-only MCL checkpoint diagnostic. Require native completion,
finite clustering weights, verified IDs and artifact provenance; a process
exit of zero alone is insufficient. Preserve failures separately from scores.

Retain mean seed-level precision, recall and F1, not pooled pair metrics.
Use 20,000 PCG64 multinomial whole-seed bootstrap replicates with seed
20261130, reset per contrast. The 14 primary F1 contrasts remain two OrthoHMM
methods versus full OrthoFinder across seven conditions, with Bonferroni
intervals; precision and recall contrasts are exploratory. Report conditional
complete-case exclusions and failure fractions, with no zero imputation.

This additional panel does not establish generalization by itself. Independent
curated validation, ablations, other robustness experiments and the remaining
publication requirements are still necessary. The protocol, mapping algorithm,
implementation, manifests and method commands must be committed and pinned
before scientific generation/inference, with any later deviation disclosed.
