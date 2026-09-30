# Frozen OrthoHMM Search And Profile Scoring

Source-level Methods specification, 30 September 2026, for revision
`7f3a9e40dd7e79f842cc2c11fb8b548f9a802806`. This does not change code, defaults,
scores or scientific admission. Historical runs retain their actual revisions;
run-specific source/backend records remain authoritative.

The [numerical readback](frozen_scoring_readback_20260930.json) hashes all nine
examined frozen source files. Tiny Python fixtures confirm transition costs,
ACDE self score 24, approximate E of 0.2623001450991833 for 1,000 target
residues, nonpositive-score handling, unknown-query emissions, strict half-gap
column removal and the MSA pseudocount formula. Three focused tests pass.
Selected frozen Python definitions run without JIT decorators; this is not
native backend, full-run, sensitivity or calibration validation. Reproduce
from the repository root with a fresh output path:

```bash
python -B -m benchmark_tools.verify_frozen_scoring --output /tmp/frozen-scoring-readback.json
python -B -m pytest -q tests/unit/test_verify_frozen_scoring.py
```

## Initial Search

High sensitivity uses BLOSUM62, exact-alphabet 4-mers, a fixed cap of 100
candidate targets per query per target species, multipass graph inference,
cluster-profile expansion and Leiden seed 4. Standard instead uses 5-mers,
an adaptive cap, no multipass/profile expansion and seed 0. Ordered species
pairs include self comparisons. The production search call requires four total
prefilter hits and one diagonal-bin hit, with diagonal-bin width 10, and
reranks candidates before scoring. These gates can omit valid homologs.

At query position i, match emissions are row q_i of integer BLOSUM62; unknown
query residues have an all-zero row. Insert emissions are uniformly -1;
the background-frequency argument does not alter these initial emissions.
Emissions are int8; transition costs and dynamic-programming scores are int32.

```text
          MM   MI   MD   IM   II   DM   DD
t = [      0, -12, -12,  -1,  -3,  -1,  -3 ]

M[i,j] = eM[i,t_j] + max(0, M[i-1,j-1] + tMM,
                        I[i-1,j-1] + tIM, D[i-1,j-1] + tDM)
I[i,j] = eI[t_j] + max(M[i,j-1] + tMI, I[i,j-1] + tII)
D[i,j] = max(M[i-1,j] + tMD, D[i-1,j] + tDD)
S = max(0, max over visited positions of M[i,j])
```

The closing cost is `max(gap_extend + 1, -1) = -1`. These are additive costs,
not normalized transition probabilities. Insertions also pay an emission cost;
deletions do not. Gap costs are not learned separately at each position.
The reference kernel uses an unreachable sentinel of -1,000,000 and two rows
per state. Unknown target residues have match emission 0 and insert emission
-1. Its half-width-64 band is centered at `floor(i * target_length /
profile_length)`; maximum length at most 50, or a nonpositive width, requests
the full matrix. Banding is not a proof of equal sensitivity to full alignment.

This is local maximum-path scoring with match/insert/delete states, not a
Forward sum over paths or an orthology posterior. The implementation alone
does not establish exact Plan7/HMMER/phmmer equivalence, despite stronger
wording in some comments. C/OpenMP/AVX2 and Numba paths exist; optional GPU
routing must not be attributed to a run without runtime evidence. This
specification does not certify all routing or cross-backend equivalence.

## Significance And Normalization

For positive raw S, the frozen implementation calculates:

```text
E = 0.134 * query_length * target_species_total_residues * exp(-0.3176 * S)
stored_search_score = S / sqrt(query_length * target_length)
```

Nonpositive S receives E = 1e10. The gate is strict `E < 1e-4`. A zero
normalization denominator is replaced by one. E is calculated before initial
length normalization; downstream helpers do not apply that same correction
again. Later graph transformations are separate from this stored score.

The significance constants come from a static matrix-specific table. Its
BLOSUM62 comment refers to gap costs 11/1, whereas the scorer uses the costs
above, including insert emissions and closing costs. This readback establishes
no empirical calibration for the actual banded score distribution. E is an
approximate filter, not demonstrated HMMER calibration, assignment confidence
or a guaranteed false-positive rate. Identical nominal cutoffs across engines
do not establish matched sensitivity or calibration. Changing this formula
would require prospective evaluation, not retrospective replacement of scores.

## Cluster Profiles

The production high-sensitivity call builds profiles for clusters of 3 through
200 genes using the built-in center-star alignment. It does not enable the
optional weakest-member jackknife or a multi-species eligibility restriction.
Failed builds are skipped. Unknown residues, gaps and dots have code 20;
columns are retained only if their gap/unknown fraction is strictly below 0.5.

For substitution matrix B, stored background q_b, observed residue count N
and unweighted empirical frequencies f_a in a retained column:

```text
R[a,b] = 2^(B[a,b]/2) * q_b / sum_c(2^(B[a,c]/2) * q_c)
r_b = sum_a(f_a * R[a,b])
p_b = (N * f_b + 5 * r_b) / (N + 5), then normalize p
eM[b] = clip(round(2 * log2(p_b / q_b)), -128, 127)
```

The zero-probability/background fallback emission is -10. Consensus is
`argmax(p)`, not necessarily raw majority. Insert emissions and transitions
remain uniform; the MSA learns position-specific match emissions, not
position-specific transitions.

Genes are prefiltered against profile consensus using 4-mers, four total hits,
one diagonal-bin hit, no reduced alphabet and at most 20 candidate profiles
per gene. The scoring half-width is 64. Profile E uses profile length and
total residues in the global gene database, not one species. Significant
profile hits retain raw S, without initial-search length normalization.

The default acceptance threshold is the minimum in-sample source-member
score. An outside gene must meet it and the E gate. The best eligible profile
per gene is selected, with the first retained row breaking a raw-profile tie.
A positive initial sequence hit to a source member must support the new edge
in either direction. Forward hits take precedence over duplicate reverse
relations. The highest remaining anchor score is used, with source-member
order breaking anchor ties; undirected deduplication follows. Raw profile
score is an acceptance/ranking statistic here, not the new edge weight.
In-sample thresholds are not calibrated assignment probabilities.

## Primary Sources And Scope

All source links identify the full frozen revision:

- [Configuration](https://github.com/JLSteenwyk/orthohmm/blob/7f3a9e40dd7e79f842cc2c11fb8b548f9a802806/orthohmm/accuracy.py) and [production call](https://github.com/JLSteenwyk/orthohmm/blob/7f3a9e40dd7e79f842cc2c11fb8b548f9a802806/orthohmm/orthohmm.py).
- [Profiles](https://github.com/JLSteenwyk/orthohmm/blob/7f3a9e40dd7e79f842cc2c11fb8b548f9a802806/orthohmm/search/profile.py), [recurrence](https://github.com/JLSteenwyk/orthohmm/blob/7f3a9e40dd7e79f842cc2c11fb8b548f9a802806/orthohmm/search/viterbi.py) and [search/normalization](https://github.com/JLSteenwyk/orthohmm/blob/7f3a9e40dd7e79f842cc2c11fb8b548f9a802806/orthohmm/search/engine.py).
- [E-value calculation](https://github.com/JLSteenwyk/orthohmm/blob/7f3a9e40dd7e79f842cc2c11fb8b548f9a802806/orthohmm/search/evalue.py) and [matrix constants](https://github.com/JLSteenwyk/orthohmm/blob/7f3a9e40dd7e79f842cc2c11fb8b548f9a802806/orthohmm/search/matrices.py).
- [MSA emissions](https://github.com/JLSteenwyk/orthohmm/blob/7f3a9e40dd7e79f842cc2c11fb8b548f9a802806/orthohmm/search/msa_profile.py) and [strict expansion](https://github.com/JLSteenwyk/orthohmm/blob/7f3a9e40dd7e79f842cc2c11fb8b548f9a802806/orthohmm/search/profile_expansion.py).

This is a numerical description, not biological validation, sensitivity
matching, uncertainty estimation or release admission. Older review archives
retain their historical manuscript bytes and do not automatically include it.
