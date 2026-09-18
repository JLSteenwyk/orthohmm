# Frozen QfO Factorial Statistics Implementation

`bootstrap_qfo_factorial.py` implements the previously frozen protocol, including
its pre-scoring native-pair clarification. Required protocol SHA256:
`f8946e12cefcf84abbee0fb9492f240c05508e045efe00a3006304d34c1fd115`.

No empirical factorial interval report has been generated. Complete admitted
eight-cell raw-count assembly remains pending on inference/scoring completion.
Passing statistical tests does not establish experiment completion or accuracy.

## Statistic and Contrasts

The input must contain all eight P/C/R cells and all18 disjoint represented
SwissTrees families, with matching reference membership/truth totals and the
frozen reference checksum. For each family, raw counts are converted by the
native raw/2+1 rule. Each replicate averages family precision and recall, then
computes their harmonic F1; it does not average family F1 or pool gene pairs.

The CLI fixes100,000 PCG64 multinomial draws at seed20260922. One common draw
matrix is used across all cells. Twelve signed on-minus-off simple effects
and two C-by-R differences of differences yield42 metric endpoints. The
adjusted percentile bounds use .05/84 and1-.05/84, with linear quantiles.
Family-level differences and positive/tied/negative counts remain descriptive;
positive interaction is not a method win. Nominal intervals are also retained.

## Verification and Limits

Fifteen new tests establish contrast orientation/inventory; independently
repeat families, recompute from raw confusion counts and reproduce all42
intervals; check exact-zero identical-cell behavior; and reject incomplete,
misordered, overlapping, mutated or inconsistent inputs. Combined with the
existing SwissTrees stage/comparator tests,36 tests pass. Synthetic fixtures
are tests only and are not written as empirical benchmark results.

The eventual CLI requires a supplied counts SHA and checks all input records
from the upstream count audit before and after analysis. Complete raw relation
identity and native-score agreement must be established by that upstream
audit; equal truth totals alone would not establish identical reference pairs.
The model remains conditional on18 development-exposed families and does not
account for method selection, shared evolutionary history, or other QfO
endpoints. No independent validation or general-superiority claim follows.

```bash
python benchmark_tools/bootstrap_qfo_factorial.py \
  --counts VERIFIED_COUNTS.json --counts-sha256 VERIFIED_SHA256 \
  --protocol benchmark_tools/results/QFO_FACTORIAL_PROTOCOL_20260917.md \
  --output NEW_INTERVALS.json --markdown NEW_INTERVALS.md
```
