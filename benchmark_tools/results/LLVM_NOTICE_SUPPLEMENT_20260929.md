# Relocatable LLVM Notice Supplement

The existing source-notice exporter now supports the recipe-bound LLVM source
inventory through an explicit `--llvm` mode. Its original four-file igraph
component selection is unchanged. The new mode exports only the top-level
LLVM, lld and compiler-rt `LICENSE.TXT` files, 46,987 bytes in total.

The exporter reads these bytes from the verified source archive, not the
inventory's copied text. It checks the recipe receipt and source URL/hash,
the notice member sizes/hashes, regular-file type, uniqueness and completeness.
Watched source identities are rechecked before finalizing the index. Failed
exports cannot overwrite an existing directory or finalize an incomplete index.

[Actual export](llvm_notice_export_20260929.json) and
[relocated verification](llvm_notice_relocation_20260929.json) both pass.
The identical 3,035-byte `SOURCE_NOTICE_INDEX.json` has SHA-256
`7614880b1fcc552d800b6ee4f98e72ceffeb1fc9968ff6b4ec31ce2116b9b01c`.
The existing offline verifier uses that externally retained digest and the
exported files; it does not need original archive paths to remain accessible.

```sh
python -B -m benchmark_tools.export_bundled_source_notices --llvm \
  --inventory benchmark_tools/results/llvm_source_notices_20260929.json \
  --output /fresh/llvm-notices --receipt /fresh/llvm-notices-receipt.json
```

Local exports are under `benchmarks/work/llvm_notice_supplement_20260929`
and `benchmarks/work/llvm_notice_relocated_20260929`. These supplement the
existing wheel notices without altering historical notice packages.

All 57 focused LLVM/source/wheel notice tests pass, including synthetic
relocation after original inputs are deleted, corrupt evidence, missing or
duplicate entries, symlinks, payload changes and extra export files.

This is selected notice packaging, not complete attribution of LLVM's
third-party components, proof that all three components are in the wheel,
source-to-binary equivalence or redistribution clearance. No source code,
source archive, binary or dataset is exported. Frozen inference settings and
runtimes remain unchanged; timing-scope approval and publication gates remain
open.
