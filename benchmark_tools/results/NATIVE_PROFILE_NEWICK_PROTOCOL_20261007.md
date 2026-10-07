# Profile Localization Saved-Tree Readback

The selected trace already exists with SHA256
`0d09567de8828b46312aa280f0bdabe6bbfbfaa01f591cbbca9190241bdb5ae0`.
Its source was committed/pushed as `247033cd` before execution. Preserve its
report, source and actual execution receipt; do not repeat the trace.

Commit/push the new `readback_allocated_native_qfo_profile_newick.py` and
focused tests before one selected execution. Require a fresh destination
`native_qfo_profile_newick_readback_20261007_v1.json` under results. Use the
existing scientific Python3.10 environment with injection variables cleared,
bytecode disabled and numerical threads limited to one. Freeze the new reader
after selected execution; any failure remains explicit, with no automatic
rerun or alteration of bound sources.

Use unchanged `readback_native_qfo_swiss_reconciliation.checked_tree` to check
checkpoint identity, raw/rooted/reconciled digests, unique checkpoint and tree
leaf universes, and identical rooted/reconciled clade topology. Compute LCAs
using Biopython/Newick, not the primary node-table LCA kernel. A separate CSV
stream of admitted root-HOG output collects complete members of selected
families, verifies root membership of both endpoints, and checks candidate
membership hashes/sizes against the trace and saved trees.

Check every changed pair and the complete trace transition inventory. For
same-family endpoints, verify the annotated Newick LCA ID/pair-event and
prediction rule; for separated endpoints, require different saved family
trees, no claimed LCA and no prediction. Verify pair/reference truth,
membership-equality flags and classification counts, plus species-tree hash
equality without assuming it. Recheck all direct bindings and newly observed
tree/checkpoint records before writing the report.

The expected selected coverage is four pairs, six family checkpoints and
70 tree leaves. Before states have speciation LCAs; after states are separated
candidate families. These are predictions to test, not fixture-derived results.
Tests use invented trees and root memberships, including same-family
duplication exclusion, corruption and fresh-output refusals. They do not
execute the selected production readback before its source commit.

Saved trees/checkpoints are newly observed current-byte records, not
originally inventoried admission evidence. Neither this reader nor the trace
establishes biological correctness, individual profile-edge causality,
calibrated node confidence or independent biological replication. There is
no new inference, alignment, reconciliation, scoring, bootstrap or timing
admission. Existing readback and frozen source files remain unchanged.
