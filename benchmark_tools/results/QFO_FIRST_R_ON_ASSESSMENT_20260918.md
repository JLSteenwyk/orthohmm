# First R-On Factorial Assessment

Cell `p0_c0_r1` completed all six native challenges. Assessment21697
completed0:0 in37:38 with8CPUs/64GiB; validator21698completed0:0 in15s.
Independent admission verifies15native tasks,48assessment records, exact
input/environment/executor bindings and output hashes. Re-running the frozen
admission command produced an identical `assessment` object; scheduler text
is time-dependent and is not treated as a byte-identical whole-report replay.

Preserved report `qfo_factorial_assessment_p0_c0_r1_20260918.json`, SHA-256
`d86179980eed61c0dd0eaac5a862498bd9280dd6489ac2b461844c878673a37d`.
This is original-release-limited QfO evidence, not a corrected-input result.

| Endpoint | Matched R-off: p0_c0_r0 | R-on: p0_c0_r1 |
| --- | ---: | ---: |
| GO mean Schlicker | 0.472537 | 0.490148 |
| EC mean Schlicker | 0.931175 | 0.968251 |
| VGNC harmonic TPR/PPV | 0.668209 | 0.896931 |
| SwissTrees harmonic TPR/PPV | 0.677567 | 0.781539 |
| TreeFam-A harmonic TPR/PPV | 0.576755 | 0.562210 |
| FAS mean | 0.775253 | 0.783598 |
| Project-defined secondary six-metric mean | 0.683583 | 0.747113 |

R-off values come from `qfo_factorial_assessment_p0_c0_r0_20260917.json`.
All displayed endpoints are retained; the TreeFam-A decrease is not omitted.
P-off still includes the initial HMM search. Native inferred pairs replace
group-clique pairs in R-on. Both cells have profile expansion and candidate
expansion disabled; this is not the satellite_v2 production comparison.

For R-on, SwissTrees native recall/precision are0.66492905/0.94774636;
TreeFam-A0.39794700/0.95740331; VGNC0.8133617448/0.9996405464.
The retained prediction count is4,950,789reference-mapped native pairs.
No interim significance claim, tuning, cell selection or bootstrap is made.
Wait for all eight admitted cells before the prespecified paired SwissTrees
factorial uncertainty analysis. TreeFam source-family recovery remains open.

Runtime above is scoring on the shared analysis host, not matched inference
timing. Native FAS sampling and endpoint-specific uncertainty semantics are
unchanged. Native `stderr` fields are not paired method-difference intervals;
in particular, GO/EC fields carry the separately audited confidence-halfwidth
semantics. The six-metric mean is not official QfO F1 or a superiority test.
