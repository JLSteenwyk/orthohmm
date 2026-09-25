# Seven-Protein Search Trace

Job **22169 completed 0:0 in 4:24** from clean frozen source
`e1e96f2d058f4cdaa5875119d0e25bd1e07a8226`. The read-only scan verified all
401,908,030 rows and 36,360,701,687 bytes against the recovered search table
SHA256 `2d395e39975cbb73c0bd33f016836aa0621e4101b378d809f2579ac37a42f117`.
Elapsed time is descriptive shared-host accounting, not comparative timing.

The [report](qfo_native_residue_hit_trace_20260925.json) and
[complete incident HSP subset](qfo_native_residue_hit_trace_20260925.m8)
are verbatim copies of the execution outputs. Report SHA256 is
`a51bf545aa924286e08a205416536596cc14895f4763ba19b45cb39aea538cb4`;
the 1,360-byte subset SHA256 is
`5229ddfa487429ed4058e88bccdb83f43cba49e9f5d49f0a253bb59ec764f99a`.

| Proteins | Outgoing partners per protein | Incoming partners per protein |
| --- | ---: | ---: |
| P58865, P58866 | 2 | 2 |
| Q8TN68, Q8TS72, Q8TTA5 | 3 | 3 |
| Q8TS73, Q8TTA9 | 2 | 2 |

Counts include self hits. Each row above is a fully reciprocal component
in this selected search graph; there are no hits in either direction between
the seven proteins and any other input protein. The 17 retained HSPs are
17 distinct directed pairs: seven self pairs and ten non-self directed pairs
(five unordered non-self pairs). Every protein has a self hit.

Independently reparsed the retained subset to check HSP counts, distinct pairs,
self hits, reciprocity, and that every endpoint belongs to the seven targets.
Matched the target identities to the earlier
[reference-exposure report](qfo_native_residue_reference_exposure_20260925.json):
all belong to `UP000002487_188937`, so there are zero observed cross-species
incident search pairs. Rechecked subset/source hashes and executor cleanliness.

This describes the retained legacy search, including its seven O deletions.
It does not establish final OrthoMCL membership, ortholog pairs, score effects,
or equivalence to a residue-preserving counterfactual. In particular, absence
of observed cross-species hits cannot exclude hits that a different residue
representation might recover. The reference-exposure and final-group analyses
remain separate from this search trace. The full publication goal remains open.
