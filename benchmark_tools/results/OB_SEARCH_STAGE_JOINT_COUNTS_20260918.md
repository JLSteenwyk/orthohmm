# Direct Search Evidence And Grouping Outcomes

This exploratory analysis summarizes all70 admitted OrthoBench family traces,
not selected examples. The40733 pair memberships include low-certainty and
within-species pairs; they are not official weighted benchmark denominators.
Overlapping reference memberships are preserved. No method is retuned.

| Retained direct-hit directions | Pair memberships | Together, profile-off refined | Together, profile-on refined | Together, candidates | Together, root HOGs |
| --- | ---: | ---: | ---: | ---: | ---: |
| None | 24333 | 6010 | 6252 | 9682 | 8465 |
| One | 1321 | 496 | 550 | 782 | 648 |
| Both | 15079 | 13015 | 13395 | 14005 | 13781 |

At the root-HOG stage,15868 no-hit pairs,673 one-direction pairs and1298
two-direction pairs remain separated. Thus1971 separated reference-pair
memberships have retained direct search evidence. Conversely,8465 no-hit
pairs are grouped together. A missing direct hit is neither necessary nor
sufficient for final separation: indirect graph paths and later grouping
operations can connect pairs without a direct retained hit.

Candidate-to-root processing separates1575 previously grouped pair
memberships:1217 with no direct hit,134 with one direction and224 with both.
It gains no co-membership pairs, consistent with root HOGs remaining within
candidate groups. This transition includes tree inference, reconciliation
and constraints; it does not identify which operation caused a split or
whether the reference relationship represents the correct biological history.

The profile-off and profile-on refined columns are matched branches, not
consecutive stages. Cached direct hits are not necessarily selected RBNH
edges, and absent hits do not distinguish prefilter rejection from scoring
rejection. This analysis does not reconstruct missing rejected-hit logs or
attribute causal search/edge failures. It supplies neither a new confidence
interval nor independent validation.

## Verification

`summarize_ob_search_stage_trace.py` consumes the admitted trace report with
SHA-256 `bda00fb593b357bc8f07e43544feae598150ae891ed82c988fb339a92349de0c`
and its pinned pair table. Every family must contain exactly its complete
unordered pair universe; duplicate/foreign pairs, malformed booleans and
nonfinite/nonpositive scores are rejected. All search and stage marginals
must reproduce the admitted report. Source hashes are checked before/after.
Eleven tests cover these checks and grouping despite absent direct hits.

`ob_search_stage_joint_counts_20260918.json` retains every family's joint
counts plus aggregate membership counts, both across all pairs and separated
into cross-species/within-species scopes. Absent counter keys mean zero.
Only admitted trace identities are rechecked here; the expensive original
native partition and hit-cache validation is not rerun by this summarizer.
