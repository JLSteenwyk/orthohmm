# Frozen YGOB Interpretation

First held-out accuracy inspection followed successful native-file, command,
conversion, independent reference-reconstruction and overlap admission checks.
The retained pair-level TP/FP/FN counts agree exactly with direct enumeration
for all four methods. No thresholds, groups, exclusions, or endpoints were
selected using these outcomes.

The primary satellite_v2 minus full OrthoFinder F1 difference is -0.084870
percentage points. Its adjusted interval [-0.622528, 0.445222] includes zero.
This is neither superiority nor a formal equivalence/noninferiority result;
no equivalence margin was prespecified. Higher precision and lower recall
are supported separately, not an unqualified accuracy advantage.

High sensitivity performs worse than full OrthoFinder on F1, precision and
recall, with all adjusted intervals excluding zero. The sequence-only
OrthoFinder checkpoint is diagnostic, not an independently timed method.
The difference between OrthoHMM modes cannot isolate reconciliation alone:
satellite_v2 also changes candidate-family expansion.

The target is curated homolog-group co-membership, including within-species
pairs, rather than resolved ortholog pairs. Every input/reference gene is
represented, so missing-prediction coverage does not explain these aggregate
differences. Thirteen input-only genes remain projected out as specified.

This is novel-taxon transfer, not family-disjoint validation. Development
homology hits occur in 85.98% of proteins and 67.82% of pillars. Curation shares
sequence and annotation ancestry; absence of label input to the reviewed
inference workflows does not eliminate this dependence. New method tuning
based on these results needs a new independent confirmation dataset.

Timing remains shared-machine and not controlled comparative efficiency.
Dataset redistribution permissions remain unresolved. Completing this
evaluation does not complete the publication objective.

See [results](YGOB_FROZEN_RESULTS_20260916.md) and
[machine-readable summary](ygob_frozen_results_20260916.json). The summary
links the full per-pillar results by path and SHA256; those larger generated
results remain preserved locally for reproducibility and archive assembly.
