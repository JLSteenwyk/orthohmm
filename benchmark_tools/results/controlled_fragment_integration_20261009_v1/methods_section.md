### Controlled Fragment Observation Test

One prospectively frozen observation condition reuses the ten variable-length
baseline histories (seeds 20261101--20261110), not new biological histories.
Within each dataset, SHA256 of UTF-8 `seed:condition:gene_id` ranks all gene IDs;
the first floor(N/5), with ID tie-breaking, retain a centered floor(3L/5)
substring starting at floor((L-floor(3L/5))/2). Identifiers, species ownership,
unselected sequences and the byte-exact parent orthology truth stay unchanged.
Independent readback verified all 8,200 genes and 1,636 truncations. This is
controlled synthetic truncation, not natural fragment annotation or a new
independent generalization set. The simulator does not add within-family indels
or validated domain architectures in this condition.

All 30 fresh native runs completed: high-sensitivity OrthoHMM (HS), satellite_v2
phylogenetic OrthoHMM (PHY), and full OrthoFinder 3.1.5 (OF), sequentially for
each seed. The OF sequence-only MCL checkpoint is an additional diagnostic
output, not a separately rerun final sequence-only workflow or timing estimate.
Forty retained baseline outcomes are reused without rerunning controls. Settings
are frozen except input/output/metrics paths. The OrthoHMM CPU budget is four
with four threads per worker; actual metadata reports one search worker and
four total search threads, not four simultaneous workers. The sequential job
received 16 CPUs and 16 GiB on the shared Threadripper. The prospective runtime
amendment records 72 full-inventory differences, while eleven required
scientific dependency versions and the full OF inventory remain unchanged.
This is not proof of historical OrthoHMM output equivalence.

Micro cross-species pair F1, precision, recall and pair-endpoint coverage are
calculated per seed; displayed means are arithmetic means across the ten
seeds, not pooled-pair F1 or a harmonic mean of mean precision and recall.
Coverage is the fraction of all input genes present in at least one predicted
cross-species pair, not recall. Zero/one/two truncated-endpoint strata use the
same prospective flags in both arms and include inter-origin-family false
positives. Their integer TP/FP/FN counts partition each whole-dataset count.

The fixed 15-endpoint comparison family comprises fragment-minus-baseline for
HS, PHY and OF plus fragment HS-minus-OF and PHY-minus-OF, each for F1,
precision and recall. Whole seeds, not dependent gene pairs, are resampled
20,000 times with PCG64 seed 20261011 and linear quantiles. Nominal intervals
use .025/.975 and fixed-15 Bonferroni intervals use .05/30 and 1-.05/30.
These are conditional percentile approximations with only ten seed units;
extreme-tail coverage is not guaranteed. All planned seed lists are complete.
No checkpoint comparison interval or subgroup interval is added after seeing
results. [Frozen protocol](CONTROLLED_FRAGMENT_OBSERVATION_PROTOCOL_20261009.md),
[runtime amendment](CONTROLLED_FRAGMENT_RUNTIME_AMENDMENT_20261009.md),
[executed validation](controlled_fragment_execution_20261009_v1.json).

