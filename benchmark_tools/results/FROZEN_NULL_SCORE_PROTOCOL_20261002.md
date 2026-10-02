# Frozen Search Null Score Protocol

Freeze before inspecting null-score outcomes. This is a supplementary synthetic
diagnostic of the documented approximate significance filter, not new benchmark
accuracy, independent biological validation, parameter optimization or timing.
The scientific configuration, original scores and all defaults stay unchanged.

## Design

Use scientific revision `7f3a9e40dd7e79f842cc2c11fb8b548f9a802806`, independently
retained scalar C kernel SHA256
`f21d902852d727d6270b65ba900db52e7216cd874ee2d2b89fa72b6f16dcd5c0`,
BLOSUM62, existing match/insert/delete transitions and static lambda=0.3176,
K=0.134. Verify exact Git source blobs, native binary and imported origins;
require the existing NumPy 2.2.6/Numba 0.65.0/llvmlite 0.47.0/SciPy 1.15.3
runtime. Do not upgrade shared packages or silently fall back to another kernel.
Use two scalar-C threads and two JIT threads. No performance claims follow.

Generate query and target residues independently, not by mutation of a shared
ancestor. Each pair constitutes its own one-target database: query length,
target length and database total residues are equal. Packed execution batches
do not combine their databases. Retain 1,000 pairs for each of ten seed indices
0..9 per cell, with PCG64 SeedSequence `[20261002, seed, regime_index, length]`.
Query draws precede target draws. Retain sequence digests and every raw score.

Lengths: 50, 150 and 400 residues. Composition regimes, in fixed order:

1. Normalized frozen BLOSUM62 residue background.
2. Uniform probability over all 20 canonical residues.
3. Half the normalized background plus probability 0.5 at glutamine, index 13
   in the frozen `ACDEFGHIKLMNPQRSTVWY` alphabet.

All nine cells have 10,000 independent pairs. For each exact pair, score both
full-matrix mode (width 0) and production width 64; width 64 uses the frozen
full-matrix fallback at length 50. The study has 90,000 independent pairs,
180,000 native score evaluations and 90 tail endpoints. The two band outcomes
for a pair are dependent; do not double their sample size. Use only scalar C
for this diagnostic, not an assertion about historical run routing. Require
banded scores not to exceed full scores and verify actual frozen E-value
decisions against exact integer-score boundaries.

## Statistics

Preserve the production strict E < 1e-4 gate. Additionally describe E cutoffs
1, 0.1, 0.01 and 0.001, without promoting them to new defaults. E is computed
from raw integer S before any length normalization:
`0.134 * length**2 * exp(-0.3176 * S)`; nonpositive scores get E=1e10.

For each cutoff, calculate the smallest positive integer score satisfying
the strict inequality and report its actual approximate E. The Poisson-model
tail reference is `1 - exp(-E_at_that_score)`, not the nominal cutoff itself.
This distinction preserves score discreteness and avoids treating E as a
probability. The [NCBI statistical tutorial](https://blast.ncbi.nlm.nih.gov/tutorial/Altschul-1.html)
explains that E is an expected count and relates it to a Poisson tail. Applying
that reference here is a diagnostic hypothesis, not a theorem for this scorer,
its gaps, band, compositions or finite lengths.

Pool the independent pairs within each composition/length cell, preserving
all per-seed scores. Report successes, denominator, observed tail fraction,
observed/model-reference ratio and exact Clopper-Pearson intervals using
[SciPy's documented method](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats._result_classes.BinomTestResult.proportion_ci.html).
Report nominal 95% intervals and Bonferroni intervals across all 90 endpoints,
with confidence level `1 - 0.05/90`. These describe Monte Carlo uncertainty
conditional on the synthetic generator, not benchmark/reference-family
uncertainty. No additional p-value, tail-fit, coefficient estimation, cutoff
selection or alternative-model search. Zero observed hits is not zero error
probability. Ten thousand trials may be insufficient at the production cutoff.

Report paired band score changes and lost threshold decisions descriptively,
retaining every cell and seed. Do not fit a new lambda/K, extrapolate lower
tails, claim a global family-wise error guarantee or tune on the null results.

## Execution And Interpretation

Commit this protocol, producer and unit tests before scoring the panel. Run
in a fresh interpreter with an explicit frozen checkout and fresh output:

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2 \
python -B -m benchmark_tools.probe_frozen_null_scores \
  --root benchmarks/work/publication_method_native_v2 \
  --protocol benchmark_tools/results/FROZEN_NULL_SCORE_PROTOCOL_20261002.md \
  --output benchmarks/work/frozen_null_score_panel_20261002
```

Recheck source/protocol/native identities at each completed cell. Retain all
partial/failed attempts, and never overwrite or silently retry. No k-mer
prefilter, learned cluster profile, graph, reconciliation, benchmark scoring,
host workload probe or controlled resource measurement is executed. Forced
scoring isolates the numerical filter; it does not measure production search
false positives after candidate selection. Strong composition is deliberately
synthetic, not an estimate of real-proteome prevalence or empirical evolution.
Whatever the outcome, calibration and generalization remain unestablished
outside this finite diagnostic. Further scientific changes require prospective
evaluation and cannot replace any historical scores or admissions.
