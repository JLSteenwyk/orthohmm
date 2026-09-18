# Corrected Sequence-Control Self-Hit Semantics

## Completed Conversion, Admission Pending

Job21791 completed0:0 in01:56:01 with2CPUs. Its conversion manifest is
benchmarks/results/qfo_sequence_numeric_v1/manifest.json, SHA-256
acc53be068064f84336c699b15979982b7ef8c16bdbfcf4aa1606fbe159a829b.
The values below are conversion-reported numerical counts, not yet an
independent source-equivalence admission or scientific accuracy result.

| Variant | Genes | Species | Directed hits | Self hits |
|---|---:|---:|---:|---:|
| all_hits | 984137 | 78 | 593510904 | 980829 |
| top100 | 984137 | 78 | 321164891 | 980803 |

The self-hit count differs by26. Neither panel reports one self-hit for every
input protein. Do not infer the identities or biological cause of absent
self-hits from these totals. Independent source validation21792 remains live.

## Frozen Ranking And Downstream Behavior

The declared cap ranks each query/target-species partition by descending raw
score, then target gene ID (indices follow lexical IDs), retaining rows1..100.
There is no self-hit exception. Thus a self-hit can fall outside the cap,
including on exact ties when its ID sorts after100 alternative targets.
Conversely a retained self-hit consumes one slot, leaving at most99 non-self
hits in that same-species partition. The statement that self-hits are retained
means they are not explicitly removed from the search input/all-hit variant;
it does not guarantee that all survive a diagnostic cap.

Direct downstream checks:

- RBNH construction excludes query==target before reciprocal-best selection.
- Singleton assignment requires singleton query and non-singleton target, so
  a self-hit cannot qualify.
- Cross-cluster refinement aggregation requires different source/target
  clusters, excluding self-hits.

Seven new regression cases cover cap/tie behavior and invariance to directly
adding/removing self rows in those three operations.18 focused tests pass in
total with the existing converter tests. This does NOT show that changing the
cap or reserving self slots would leave clustering unchanged: those changes
can alter which non-self hits survive. No such change has been made.

## Source Identity And Limits

The tested accuracy.py and refinement.py bytes match publication core
7f3a9e40dd7e79f842cc2c11fb8b548f9a802806:

- accuracy.py:1a35944ab7fea859143f599b1262aa111272787f6f735b9acc37d299427f2ad6
- refinement.py:991f1eb6a5f73d0442529ed19095a34b7c6ba8bff8dfe43a1127e24ec73fb26d
- converter helper:b451d89a53102128dee4e65b64b7b970a4a048d51bec1bc6c4aabf762597f70c,
  matching the helper identified by the completed conversion manifest.

Exact source/checkpoint equality remains the responsibility of21792. This
review establishes code semantics and a numerical observation, not which
specific26 self rows were excluded or why their search scores ranked as they
did. Do not change the prespecified diagnostic or claim matched HMM sensitivity.
Graph memory review21798 still waits for independent admission.
