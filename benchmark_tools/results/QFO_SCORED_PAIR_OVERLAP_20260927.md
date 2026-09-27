# Corrected GO/EC Scored-Pair Overlap

Descriptive analysis of the first three methods in the frozen corrected
comparison v7 manifest (SHA256
`042aa221d01554114a1b0b413ca3e2ff56dad5f8c06362c69a1bf49f2de642fc`).
The manifest, admissions and endpoint records were checked before/after
analysis. Each raw count matches its admitted assessed-relation count, and
each rounded mean matches its admitted mean within 5.051e-7. This uses
corrected runs, not the older GO/EC arithmetic audit's historical predictions.

| Endpoint | Left | Right | Left pairs | Right pairs | Shared pairs | Shared scores differing |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| GO | OrthoHMM high | OrthoFinder full | 145,619 | 163,557 | 78,352 | 0 |
| GO | OrthoHMM phylogeny | OrthoFinder full | 84,211 | 163,557 | 75,894 | 0 |
| GO | OrthoHMM phylogeny | OrthoHMM high | 84,211 | 145,619 | 74,994 | 0 |
| EC | OrthoHMM high | OrthoFinder full | 185,664 | 175,361 | 120,158 | 0 |
| EC | OrthoHMM phylogeny | OrthoFinder full | 117,460 | 175,361 | 101,712 | 0 |
| EC | OrthoHMM phylogeny | OrthoHMM high | 117,460 | 185,664 | 110,604 | 0 |

All shared-pair scores are identical at their retained six-decimal precision.
Aggregate differences therefore arise from differing scored-pair membership
and denominators at this precision, not differing scores for shared pairs.
This does not prove identity of unavailable full-precision values, underlying
annotations or native scorer implementation across runs.

The [summary](qfo_scored_pair_overlap_20260927/summary.json) and six adjacent
reports retain exact integer-millionth score sums, pair counts, rounded means,
input hashes and an additive decomposition of each original mean difference:
shared-score sums divided by each method's original count, plus left-only
contribution minus right-only contribution. Shared conditional means are
diagnostic only: replacing the native endpoints with an intersection-only
mean would change the estimand and erase these observed differences.

## Reproduction And Limits

Use the ordered raw-file inputs recorded in each report:

```bash
python -m benchmark_tools.compare_qfo_scored_pairs \
  --left LEFT_RAW.txt.gz --right RIGHT_RAW.txt.gz --metric GO \
  --output NEW_REPORT.json
python -m pytest -q tests/unit/test_compare_qfo_scored_pairs.py
```

Fourteen tests cover parsing, canonical pair direction, duplicate/invalid
records, exact serialized sums, disjoint/identical sets, distinct denominators
and output overwrite refusal. No native inference or scoring was rerun.
These are eligible scored pairs, not prediction coverage over the entire
proteome or independent family units. No confidence intervals, causal effects,
new default selection or biological superiority are established. The lack of
appropriate independent units remains an uncertainty-analysis requirement.
