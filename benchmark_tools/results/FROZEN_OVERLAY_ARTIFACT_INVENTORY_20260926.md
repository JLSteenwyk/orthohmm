# Frozen Overlay Artifact Inventory

Reused the tested wheel/notice readers on the exact artifacts in the September
26 setup-overlay installation, rather than assuming that the older development
wheelhouse inventory covered the new compiled runtime. The installation report
is pinned by SHA-256
`00f6a3457d32ce860e4f6afc6d924e53a09a8663ee99f5d150bbe2f526a578eb`
in the frozen-source installation evidence.

## Results

- All 11 wheel hashes, distribution/version names and embedded metadata match
  the installation report. There are 79 notice candidates, 47 native-library
  members and zero unresolved `License-File` declarations.
- All ten third-party wheels are byte-identical to the earlier development
  wheelhouse. The rebuilt OrthoHMM wheel is different, as expected for the
  frozen-source/setup-overlay artifact. Its notice candidates, declarations
  and native-member names are unchanged. Its metadata headers are unchanged;
  the metadata body differs, so complete metadata identity is not asserted.
- All 42 entries of the new OrthoHMM wheel are covered by its `RECORD`, with
  exact recorded hashes and sizes. The embedded project license matches
  `source/LICENSE.md` from the frozen staging tree.
- Each of the three project native libraries declares `libgomp.so.1` and
  `libc.so.6`. None contains the checked `cudaMalloc`, `cudaFree` or
  `__cudaRegisterFatBinary` symbol witnesses. This bounded symbol check is not
  proof of absence of every possible static dependency or GPU-related code.
- All 30 relevant notice, wheel-content and frozen-installation tests pass.

## Evidence

[Exact notice inventory](publication_frozen_overlay_notices_20260926.json),
SHA-256 `6c4457cb2464b6ebb1a32790099c6e5cbc185a995a3045413e105cd0088b670f`;
[comparison with prior wheelhouse](publication_frozen_overlay_notice_comparison_20260926.json);
[project wheel RECORD/license/native-linkage readback](publication_frozen_overlay_wheel_contents_20260926.json),
SHA-256 `cd706c11d7586985d9e0fc6cf3209653b61f1099ea558af3cbacdd98f9b8c7c8`.

Notice inventory reproduction:

```bash
python -m benchmark_tools.inventory_dependency_notices \
  --install-report benchmarks/work/publication_frozen_overlay_20260926/install_report_v2.json \
  --report-sha 00f6a3457d32ce860e4f6afc6d924e53a09a8663ee99f5d150bbe2f526a578eb \
  --output benchmarks/work/frozen_overlay_notices_new.json
python -m benchmark_tools.audit_wheel_contents \
  benchmarks/work/publication_frozen_overlay_20260926/wheels/orthohmm-0.5.0-cp310-cp310-linux_x86_64.whl \
  --license benchmarks/work/publication_frozen_overlay_20260926/source/LICENSE.md
```

The second command prints the wheel audit; the retained report also binds the
installed-source evidence, source reader and project-license identities. No
package was installed, loaded or modified by this inventory. Native bytes were
read through `readelf` and `nm`, not executed. No third-party wheel or notice text
was newly redistributed. Full transitive linkage, provider notice completeness,
legal compatibility, compiler/OS obligations and final selected-archive review
remain open. No new score, scientific default or publication-readiness claim
follows from this artifact inventory.
