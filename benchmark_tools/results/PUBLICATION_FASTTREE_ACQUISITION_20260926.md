# FastTree Release Source and Binary Provenance

Acquired six files from the immutable upstream revision
`29c5e62fbcd93230ee325f9c6a17b81f00e3c72a`, the peeled `v2.2.0` tag:
`FastTree`, `FastTree.c`, `LICENSE`, `README.md`, `ChangeLog.txt` and
`index.html`. The acquisition script pins every SHA-256 and retrieves the
files over HTTPS without executing them. All bytes remain local under
`benchmarks/work/publication_fasttree_source_20260926/acquired/`.

The [receipt](publication_fasttree_acquisition_20260926.json) has SHA-256
`c553af49b95153434c6c09e4bf4dd6bf215d472a593c6f13103394785570ec68`.
All eight source/binary/script records were rechecked after acquisition. The
installed binary record also matches the earlier
[external-tool relocation](PUBLICATION_RELOCATED_PHYLOGENY_TOOLS_20260926.md)
record exactly. Its 1,496,928 bytes are identical to the upstream release
executable, SHA-256
`55a9d997813aae2208bd4c2081bfa690e0ecdba2d6c491805d8689415c43e38e`.
This resolves the previously missing upstream source/notice association for
that executable; it does not retrospectively prove every historical job used it.

## Reacquisition

From the repository root, with a new destination:

```sh
/usr/bin/python3 -S -m benchmark_tools.acquire_publication_fasttree \
  --output benchmarks/work/publication_fasttree_source_20260926/acquired \
  --installed /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/FastTree_v220/FastTree
python -m pytest -q tests/unit/test_acquire_publication_fasttree.py
```

Existing output paths are refused. Downloads have a size limit and timeout;
hash mismatches fail closed while retaining downloaded bytes and a failure
receipt. The installed executable is checked before and after acquisition
and is never changed. Seven tests cover successful acquisition, overwrite
refusal, wrong installed bytes, download hash mismatch/failure reporting,
oversize downloads, symlinks and insecure redirects. No scientific inference,
new benchmark score, timing run or DGX access occurred.

## Source Notices and Remaining Scope

The pinned [source header](https://github.com/morgannprice/fasttree/blob/29c5e62fbcd93230ee325f9c6a17b81f00e3c72a/FastTree.c)
declares GPL version 2 or later. The separate
[LICENSE](https://github.com/morgannprice/fasttree/blob/29c5e62fbcd93230ee325f9c6a17b81f00e3c72a/LICENSE)
contains GPL version 3 text. Preserve both notices rather than treating the
repository's license label as a replacement for the source header. This
inventory is not a legal opinion or a blanket archive-compliance finding.
No third-party binary/source files were added to this project's Git history.

The pinned upstream installation documentation identifies its Linux binary
as requiring AVX2. Matching that binary does not make it baseline-CPU portable.
The upstream changelog reports GCC 15.1.0 for Linux executables; exact build
inputs, compiler flags and dependency closure have not been independently
reconstructed. A new source build must be identified separately, tested and
must not silently replace the benchmark executable. Mutable documentation
links in the upstream README are not used as acquisition targets.

Hashes and HTTPS establish observed content identity, not a signed maintainer
attestation. The full source-build/environment reproduction, MAFFT extension
notice review, other runtime dependency sources and final artifact-level
redistribution review remain open. Controlled scaling and unresolved scientific
analyses are unaffected; the publication objective remains incomplete.
