# External Phylogeny Notices Exported And Verified

Executed the [fixed plan](PHYLOGENY_NOTICE_EXPORT_PLAN_20260930.md) once after
committing/pushing exporter, verifier amendment and 91 passing tests as
`f707e041`. Used system Python 3.12 with site initialization disabled; no
package installation, source acquisition, compilation or inference occurred.
Existing source archives, tool installations and old notice exports are intact.

## Retained Selection

| Artifact text | Bytes | Exported path |
| --- | ---: | --- |
| MAFFT 7.525 core notice | 1,763 | `mafft/license` |
| MAFFT extension notice, contextual only | 4,894 | `mafft/license.extensions` |
| FastTree 2.2.0 separate LICENSE | 35,149 | `fasttree/LICENSE` |
| Exact FastTree.c leading comment | 1,859 | `fasttree/source-header.txt` |

Total selected payload: **43,665 bytes**. The identical 5,489-byte
`SOURCE_NOTICE_INDEX.json` has SHA256
`ee32a6494607d83caf184ada6a2e197dfbeded05e139887ad7049035c8bdd1a9`.
The [byte-identical retained index](phylogeny_notice_index_20260930.json) contains
member identities, provenance, source-header interval `[0, 1859)` and source
identity, not the notice texts or code/binary payloads.

## Actual Checks

- [Export receipt](phylogeny_notice_export_20260930.json): four files verify;
  all externally supplied receipt/artifact pins and selected notice identities
  agree. All MAFFT archive names/types/duplicates/size bounds are checked.
- [Fresh-copy receipt](phylogeny_notice_relocation_20260930.json): the existing
  verifier reproduces the same index digest at a different root. A Python audit
  hook rejects reads of the five literal original artifact paths; seven observed
  verification opens use only the relocated index/notices. This is not an
  OS sandbox or proof against every possible path alias. Unit tests separately
  verify after deleting synthetic original input directories.
- [Independent readback](phylogeny_notice_readback_20260930.json): system Python
  with isolated mode and no exporter imports checks all five input identities,
  six helper Git bindings to `f707e041`, original tar/file bytes and both exports.
  A separate line-based comment-boundary calculation matches the exact 1,859-byte
  header and recorded interval; all four payload identities agree.
- All **91 focused notice inventory/export tests pass in 1.19s**, including
  corrupted evidence, duplicate/missing notices, unsafe archives, symlinks,
  changing watched inputs and FIFOs. Duration is not controlled timing evidence.

Local exports are `benchmarks/work/phylogeny_notice_supplement_20260930` and
`benchmarks/work/phylogeny_notice_relocated_20260930`. Raw notice texts remain
local; only scripts/tests/small identity receipts are committed.

## Reproduction And Limits

Run with fresh output/receipt paths, from the repository root:

```bash
/usr/bin/python3 -S -B -m benchmark_tools.export_phylogeny_notices \
  --mafft-build benchmark_tools/results/publication_mafft_build_20260926.json \
  --mafft-build-sha 4c4e92a29c1dc4b27e4b47f4010e9ea9648c17beb02b5a00995a59535394ac76 \
  --fasttree-acquisition benchmark_tools/results/publication_fasttree_acquisition_20260926.json \
  --fasttree-acquisition-sha c553af49b95153434c6c09e4bf4dd6bf215d472a593c6f13103394785570ec68 \
  --output /fresh/phylogeny-notices --receipt /fresh/phylogeny-notices.json
```

Verify any relocated copy with
`benchmark_tools.export_bundled_source_notices.verify(directory, index_sha256)`;
only its index and selected texts are read. The source-path provenance in the
index is historical metadata, not required input to that verification.

This completes selected-notice packaging, not all per-file attribution,
corresponding-source obligations, compiler/OS provenance, compatibility or
redistribution clearance. Optional RNA engines are not admitted by including
their contextual notice. FastTree's source declaration and separate LICENSE
are retained without conflation. No source code/binary/dataset archive was
published, and no scientific default or benchmark score changed. Controlled
timing, original TreeFam sources/remaining uncertainty, high-CPM admission,
final release/rights review and public deposition remain open.
