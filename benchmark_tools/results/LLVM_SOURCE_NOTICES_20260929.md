# Recipe-Bound LLVM Source And Notices

Following the [resolved-package inspection](LLVMDEV_PACKAGE_RECIPE_20260929.md),
downloaded the LLVM 22.1.0 source archive from the URL in its embedded recipe.
The 167,040,408-byte archive matches recipe SHA-256
`25d2e2adc4356d758405dd885fcfd6447bce82a90eb78b6b87ce0934bd077173`.
The [machine-readable receipt](llvm_source_notices_20260929.json) binds the
recipe receipt, source URL and local archive, and retains notice texts and hashes.

## Inspection

Streamed the archive without extracting archive paths, following links,
executing scripts or changing an environment. Checked the unique expected root,
duplicate and traversal names, a 500,000-member limit, a 5 GB total regular-file
limit and bounded reads for selected texts. Rechecked source and recipe identities
after scanning. Observed 184,760 entries: 168,946 regular files, 15,795 directories
and 19 symlinks. Regular members total 2,023,461,242 bytes.

Retained 39 filename-based notice candidates, including `llvm/LICENSE.TXT`,
and the upstream `AddLLVM.cmake` and `LLVMVersion.cmake` texts. The source's
LLVM license text includes Apache 2.0 with LLVM exceptions, third-party notice
guidance and legacy NCSA terms. The package metadata's short `NCSA` field is
therefore not an adequate substitute for these actual texts.

This closes acquisition of the recipe-named source archive and retention of
its filename-selected notices. It does not establish which third-party source
components reached the wheel, cover terms embedded only in source headers,
establish source-to-binary equivalence, or provide redistribution clearance.
Retain the embedded build script's two source edits when reconstructing this
recipe; the unmodified archive alone is not the complete build input.

Raw archive remains local in `benchmarks/work/llvmdev_manylinux1_20260929/`.
No inference or scoring was rerun. No runtime, default or publication claim changed.
