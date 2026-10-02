# Frozen Synthetic Null-Score Diagnostic

Status: completed supplementary diagnostic, not significance calibration,
orthology accuracy, production timing or publication readiness. Scientific
code, defaults and all biological benchmark scores remain unchanged.

## Prespecified Execution

The [protocol](FROZEN_NULL_SCORE_PROTOCOL_20261002.md) and producer were committed
as `186bd5e47ecaa9361a86064c33d7a57c70ce1ce3` before any study pair was scored.
The frozen scientific revision is `7f3a9e40dd7e79f842cc2c11fb8b548f9a802806`.
Three compositions, three lengths, ten seeds and 1,000 pairs per seed yield
90,000 independently drawn query-target pairs. Each pair has its own one-target
database. The scalar C kernel scores both full and width-64 bands: 180,000
evaluations, not 180,000 independent observations. Kernel, source, runtime and
protocol pins are retained. No k-mer filtering or biological inference ran.

The [compressed observations](frozen_null_score_observations_20261002.json.gz)
retain all raw integer scores, sequence digests and all 90 tail endpoints,
including every prespecified cutoff, not only favorable outcomes. Decompression
exactly reproduces the original 2,657,291-byte result (SHA-256
`343bc75bbceb59ec80b460f2bcf4080f5c1684e10354284ecfc9604b7fdba45a`).
The 172,690-byte gzip is synthetic-score evidence, not an executable runtime
archive or raw-proteome bundle. Absolute recorded origins are historical
provenance, not portable input locations. See the [receipt](frozen_null_score_receipt_20261002.json).

## Fixed Production Cutoff

Each row has 10,000 independent pairs. A hit means only that the forced pair
score passes the frozen approximate filter `E < 1e-4`. It is not a predicted
ortholog or an observed pipeline false positive.

| Synthetic composition | Length | Full hits | Width-64 hits | Changed scores | Width-64 lost hits |
| --- | ---: | ---: | ---: | ---: | ---: |
| Frozen BLOSUM background | 50 | 0 | 0 | 0 | 0 |
| Frozen BLOSUM background | 150 | 2 | 1 | 2,482 | 1 |
| Frozen BLOSUM background | 400 | 0 | 0 | 6,573 | 0 |
| Uniform residues | 50 | 1 | 1 | 0 | 0 |
| Uniform residues | 150 | 7 | 5 | 2,704 | 2 |
| Uniform residues | 400 | 12 | 3 | 6,616 | 9 |
| Half background plus half glutamine | 50 | 8,983 | 8,983 | 0 | 0 |
| Half background plus half glutamine | 150 | 10,000 | 10,000 | 4 | 0 |
| Half background plus half glutamine | 400 | 10,000 | 10,000 | 67 | 0 |

The glutamine regime allocates half the probability to glutamine and the
remaining half to the normalized background; its total glutamine probability
therefore exceeds 0.5. It is an intentionally strong synthetic composition,
not a measured real-proteome prevalence or validated evolutionary generator.

For width 64, the glutamine tail fractions are 89.83%, 100% and 100%; their
90-endpoint Bonferroni-adjusted exact binomial intervals are approximately
[88.7484%, 90.8453%], [99.9181%, 100%] and [99.9181%, 100%]. These are Monte Carlo
uncertainty for the specified generator, not family-level biological uncertainty.
At 50 residues the frozen short-sequence fallback makes both bands identical.

The integer boundary changes the model-tail reference: use `1 - exp(-E)`
at the first passing integer score, not the nominal cutoff or a posterior
orthology probability. This transformation is motivated by the
[NCBI statistical tutorial](https://blast.ncbi.nlm.nih.gov/tutorial/Altschul-1.html),
but applying it to this actual recurrence is a diagnostic hypothesis, not proof
of its calibration. The retained reference probabilities are about 0.00008022,
0.00007817 and 0.00008267 at lengths 50, 150 and 400, respectively.
The exact intervals follow the
[SciPy Clopper-Pearson method](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats._result_classes.BinomTestResult.proportion_ci.html).

The few or absent ordinary-background hits do not establish rare-tail
calibration. For zero hits, the adjusted upper bound is about 0.00081853,
which exceeds `1e-4`. The strong glutamine outcomes demonstrate that the fixed
significance formula is not universally calibrated under these synthetic
conditions. They do not establish real-data orthology false-positive rates,
explain any benchmark difference causally or show an improvement over OrthoFinder.
No coefficients, thresholds, bands or defaults are fitted or promoted from
these observations.

## Verification And Limits

The independent audit regenerated all 180 query/target digests, recomputed all
90 endpoint counts and intervals, and checked paired band changes/lost decisions
on all 90,000 pairs. Intervals were recomputed with beta quantiles rather than
the producer's binomial-result API; this is still SciPy, not an independent
statistical engine. Reference Python independently matched the first pair of
every seed cell in both bands: 180 sparse score checks, not all 180,000 native
scores. No native experiment was rerun.

The first audit attempt failed before reading results because the script import
path omitted the repository. Correcting only that path allowed the audit to
finish. Its failed attempt is disclosed in the receipt, not silently counted
as successful scoring. The one-off audit source/log receipts remain local with
their hashes; the compressed observations and retained-result tests provide
public numerical readback without depending on those local files.

Producer/scoring unit checks passed 28 cases before execution. After saving
observations, all 38 producer/scoring/readback cases pass in 1.12s, with zero
errors, failures or skips; the final test receipt is pinned in the metadata.
Run the
retained-result checks without loading the native kernel:

```bash
python -m pytest tests/unit/test_frozen_null_score_result.py tests/unit/test_probe_frozen_null_scores.py tests/unit/test_verify_frozen_scoring.py
```

This new supplement is not yet integrated into the dated manuscript render or
older review archives. Controlled resource measurements remain deferred; raw
data rights, other-QfO uncertainty, final reconciliation and a complete versioned
executable release remain open.
