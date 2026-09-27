# Filtered QfO Raw-Tree Cache Fixture

Implemented a cache selector/copier for the frozen canonical-order downstream
comparison. It requires an independently admitted output inventory, unchanged
tool identities, a complete schema-2 checkpoint, identical family membership
and the recomputed native tree-input hash on the target sequences/configuration.
Changed or absent families are excluded, and missing/unadmitted files fail.

Only four artifacts per eligible family are copied: checkpoint, raw gene tree,
candidate FASTA and alignment. Copies are independent files, not hard links.
Species-tree caches, rooted/reconciled trees and final predictions are not
copied. Excluding changed families also avoids stale checkpoints/raw trees for
a family that becomes a singleton in the canonical arm.

The [fixture receipt](qfo_reused_phylogeny_fixture_20260927.json) records actual
isolated installed execution on the retained 16-gene four-species fixture:

- One raw-tree checkpoint hit, zero remapped hits, fresh species-tree inference.
- Three root groups, 36 native pairs, four duplications and three speciations.
- All four scientific readers pass.
- Root-group file bytes and all native pairs/confidences exactly match the
  fresh fixture run. Original source artifacts remain unchanged.

Forty-five focused cache/readback/launcher tests pass. Cache tests cover altered
tools/files, changed membership and sequence/configuration hashes, omitted
artifacts, unsafe filenames, duplicate/incomplete manifests and overwrite
rejection. This is a small same-input fixture, not full QfO reuse admission or
evidence that the changed canonical candidates yield equivalent predictions.

Full retained inference 22329 remains running, with audit 22330 waiting on
successful termination. No canonical native run is submitted. Next implement
the full canonical launcher with both scheduler/admission gates and source/
artifact checks before copying anything. Do not modify the 731 sources pinned
for the queued readback. These cache helpers are new, unpinned files.
Historical scores remain unchanged; the DGX remains deferred.
