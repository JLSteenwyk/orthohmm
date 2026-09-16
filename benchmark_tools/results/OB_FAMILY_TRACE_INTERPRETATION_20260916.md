# Retained OrthoBench Family Traces

## Verified Scope

All70 reference families and40,733 within-family unordered pair rows are
retained. Reference membership has1,945 entries across1,944 distinct genes;
FBpp0309618 belongs to bothRefOG021 andRefOG068. Neither assignment was dropped.
The4,390,860-byte pair table remains outside Git with SHA256
`aabf71d4f79aff06b18890b49996895cd729bd6ac1f6fb0b7f14e79ab5295e4b`.

Job21313 completed0:0 in41s. All six checkpoints partition251,378 input genes.
Every pair's group membership and species classification was subsequently
checked against native files and FASTAs; complete pair inventories and fresh
full-reference scores were reproduced. Applying all8,440 logged merges
reconstructs the candidate partition;168 logged events touch reference genes.
These checks establish bookkeeping consistency, not biological correctness.

The verified snapshot is `ob_family_trace_verified_20260916.json`, SHA256
`bda00fb593b357bc8f07e43544feae598150ae891ed82c988fb339a92349de0c`.
[All-family table](OB_FAMILY_TRACE_ALL_20260916.md) includes neutral/adverse cases.

## Source-Defined Branches

Profile expansion starts from **unrefined multipass groups**. Multipass
refinement is a separate profile-off branch, not the input to profile expansion.
The first extraction's linear checkpoint ordering was corrected after frozen
source inspection. Original output and failure records remain available; the
pair memberships and source predictions were not changed.

| Comparison | Retained Pairs | Gained Pairs | Lost Pairs | Scope |
| --- | ---: | ---: | ---: | --- |
| Multipass to refined multipass | 18,298 | 1,223 | 0 | Profile-off refinement |
| Multipass to strict profiles | 18,145 | 567 | 153 | Profile expansion plus graph reclustering |
| Strict profiles to refined profiles | 18,712 | 1,485 | 0 | Profile-on refinement |
| Refined multipass versus refined profiles | 19,361 | 836 | 160 | Matched branch endpoints, not a direct step |
| Refined profiles to candidates | 20,197 | 4,272 | 0 | Candidate expansion |
| Candidates to root HOGs | 22,894 | 0 | 1,575 | Tree inference, reconciliation and constraints |

These are raw reference co-membership counts, including within-species pairs
and low-certainty assignments, not official precision/recall. Counts can overlap
across families. They do not show that every gained pair is correct or every
lost pair is an error under benchmark or phylogenetic orthology semantics.
The matched branch comparison demonstrates that profile-associated changes
are not uniformly beneficial even at the descriptive membership level.

## Frozen Illustrations

All six illustrations were chosen by the prespecified hash ranking within
feature strata, not by their outcomes. Numbers below are within-family pairs.

| RefOG | Genes | Possible Pairs | Refined Profiles | Candidates | Root HOGs |
| --- | ---: | ---: | ---: | ---: | ---: |
| 005 | 34 | 561 | 468 | 468 | 468 |
| 014 | 44 | 946 | 195 | 308 | 213 |
| 021 | 125 | 7,750 | 1,111 | 2,931 | 2,125 |
| 024 | 51 | 1,275 | 1,082 | 1,275 | 1,275 |
| 038 | 9 | 36 | 36 | 36 | 36 |
| 067 | 13 | 78 | 66 | 66 | 66 |

RefOG024 retains the candidate-stage increase in within-family membership.
RefOG014 andRefOG021 lose some candidate-stage pairs during final processing.
The other three examples have unchanged within-family counts across these
stages; this does not imply unchanged contamination or complete equivalence.
Per-group other-reference memberships, unlabelled-member counts and incident
merge indices remain in the full snapshot for downstream inspection.

## Remaining Mechanistic Work

The search cache contains accepted normalized hits, not a record of every
candidate rejection. Absence cannot distinguish prefilter loss from scoring
rejection. This admission verifies cached-result hashes and extraction
provenance, not an independent rerun of search. Added profile-edge identities
and initial RBNH edges remain untraced in this panel.

Inspect the relevant gene trees, rooting decisions and membership constraints
before attributing root-HOG splits to a particular reconciliation mechanism.
Independent duplication-history, fragment and domain evidence remains required.
These development-exposed traces are neither a new independent validation set
nor the prespecified biological application, and do not justify changing defaults.
