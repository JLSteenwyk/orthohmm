# External Phylogeny Notice Export Plan

Export selected notices from already acquired, exact artifacts; no network
download, compilation, tool execution, scientific rerun or binary redistribution.
This extends the existing notice-only packaging workflow, not legal clearance.

## Fixed Inputs And Selection

- [MAFFT core-build receipt](publication_mafft_build_20260926.json), SHA256
  `4c4e92a29c1dc4b27e4b47f4010e9ea9648c17beb02b5a00995a59535394ac76`.
- [FastTree acquisition receipt](publication_fasttree_acquisition_20260926.json),
  SHA256 `c553af49b95153434c6c09e4bf4dd6bf215d472a593c6f13103394785570ec68`.
- MAFFT 7.525 archive SHA256
  `2876f4adc1a2de4ed206bc40896763bf208bf1a02bda52f8bfdd91cf52d73e4a`:
  top-level `license` and `license.extensions` only.
- FastTree revision `29c5e62fbcd93230ee325f9c6a17b81f00e3c72a`:
  complete `LICENSE` and the exact first comment in pinned `FastTree.c`.
  Preserve the comment's byte interval; do not reinterpret its declaration
  as identical to the separate license file or include subsequent code.

The online [MAFFT license](https://mafft.cbrc.jp/alignment/software/license.txt)
was inspected, but it has a different copyright year from the retained archive
text; it must not replace the archive bytes. The
[immutable FastTree license](https://raw.githubusercontent.com/morgannprice/fasttree/29c5e62fbcd93230ee325f9c6a17b81f00e3c72a/LICENSE)
likewise supports source identification, not a blanket compatibility verdict.

## Execution And Verification

Commit/push the exporter and tested verifier amendment before the actual export.
Use a fresh local output directory and exclusive-create result receipt. Bind
receipt/source/archive hashes, source URLs and expected acquisition/build states.
Stream the MAFFT archive, checking all member names/types, duplicate names and
member/total-size bounds. Never extract archive paths. Recheck watched inputs
before finalizing the index. Failures cannot overwrite prior outputs or finalize
an incomplete export.

Copy the four notices and index to another fresh root. Verify there using the
externally retained index digest and the existing standard-library verifier;
it must not consult paths listed as original provenance. Retain payload hashes,
index identity and actual verification output. Original local acquisitions and
historical notice exports remain unchanged. Raw notices remain local; commit
scripts, tests and small identity/result receipts, not binaries or source code.

The shared verifier additionally rejects nonregular index/payload/extra entries
before attempting reads, including FIFOs that could otherwise hang validation.
All 91 focused notice inventory/export tests pass, including the new external
selection and offline relocation controls. This does not admit scientific
results, comprehensive security, transitive attribution or redistribution.
The optional RNA extension notice is context; those engines were not built
or verified by this export. Timing and all other publication gates remain open.
