# Sequence-Search Control

This exploratory control was specified after development benchmark and YGOB
outcomes were inspected. It does not change the frozen publication method and
cannot serve as independent confirmation. Never use YGOB to select settings.

## Search Specification

Start with the frozen 12-species/251,378-protein OrthoBench panel used by the
factorial. Build one DIAMOND2.1.11 database per target species and query all
proteins against each. This preserves species-specific database scope for
E-values without repeating 144 separate query processes. Search sequentially
with32 CPUs, very-sensitive mode, E-value1e-4, BLOSUM62, gap penalties11/1,
composition mode1, masking1, all targets and one HSP per query-target.
Retain both directions and self hits. Do not apply identity or coverage cuts.

Output qseqid, sseqid, qlen, slen, raw score, bit score and E-value. Convert
raw scores using the frozen HMM pipeline's geometric-mean length divisor,
exactly once. This holds the normalization formula constant, not the score
distribution: alignment transitions, compositional correction, significance
calibration and candidate heuristics differ. No fitted scale factor or
endpoint-guided threshold adjustment is permitted.

DIAMOND documents raw-score output, all-target reporting with -k0 and the
effect of target caps on search heuristics. The latter is why all targets are
retained for subsequent transparent diagnostics, rather than assuming a
reporting cap reproduces OrthoHMM's pre-scoring candidate cap.
[Official options](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options).

## Downstream Comparison

Primary control: initial-search replacement with profile expansion off,
candidate expansion off and reconciliation off, compared with frozen p0_c0_r0.
This is the HMM-free search/graph control. Preserve Leiden CPM0.1/seed4,
RBnH/singleton graph construction and production refinement. Additional
reconciliation or profile-expansion branches must be explicitly labeled;
turning profiles back on is no longer an HMM-free method.

Use all retained significant hits for the principal control. Also report a
diagnostic top100 per query/target-species subset, ordered by descending raw
score then target ID, constructed after search. This is not equivalent to
the HMM engine's prefilter cap. Report hit-set intersection, directional
coverage, queries without hits and threshold distributions before interpreting
accuracy. Equal E-values do not establish matched sensitivity.

Use the same official OrthoBench score and paired-family bootstrap (20,000
replicates, seed20260918). Treat two search variants times F1/P/R as six
exploratory endpoints with Bonferroni correction. Preserve negative findings.
Compare cost and sensitivity descriptively; do not claim matched efficiency
from concurrent shared-machine runs or historical cached search timings.

Before QfO extension, verify the adapter's ID/index/length/score semantics and
memory behavior. Do not launch inference from incomplete target searches.
QfO extension and additional controls remain required where feasible; this
OrthoBench experiment alone cannot establish general HMM superiority.
