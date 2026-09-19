# Search Decisions And Final Grouping

Joined the independently recounted search diagnostic to the retained final
root-HOG pair trace. All 40,733 pair memberships across 70 reference families
are retained, including same-species and low-certainty pairs. Family labels
and overlap are preserved; counts are not official weighted recall.

| Observed directional search state | Grouped in final root HOGs | Separated |
| --- | ---: | ---: |
| At least one accepted direction | 14,429 | 1,971 |
| Both directions excluded by prefilter | 8,386 | 15,129 |
| Both directions scored but not significant | 32 | 218 |
| One prefilter exclusion, one non-significant score | 47 | 521 |

The three no-accepted-hit categories recover the earlier trace's 8,465
grouped and 15,868 separated pairs. Most absent-direct-hit pairs were
prefilter-excluded in both directions, but 8,386 such pairs were nevertheless
grouped together. An accepted direct hit was also not sufficient for final
co-membership. Indirect paths and later processing prevent causal attribution
of final accuracy differences to a single rejection stage.

The [machine-readable summary](ob_search_grouping_join_20260918.json) retains
same-species/cross-species subdivisions and all 70 family summaries. Its
SHA-256 is `ecd37950fd17e1221382ad25e930e5f87bd91d424e30632658ec680fbdc9d1f8`.
The full joined table is retained at
`benchmarks/work/ob_search_grouping_join_20260918/pairs.tsv`, SHA-256
`beaaf9eb09409d5f9eb2198a40b09a988d48f8ec30754c07c67e1ff80b6489b8`.
An independent AWK recount of that table reproduced all eight collapsed
cells above and the 40,733 total. This checks arithmetic, not independent
biological truth or native-score computation.

The join rehashes the pinned admission, driver report, direction tables and
original family trace, checks exact per-pair historical hit presence, and
rejects duplicate family pairs. Thirteen focused tests pass. No search,
clustering or phylogeny was rerun, no prediction changed, and no hypothesis
test or confidence interval was introduced.
