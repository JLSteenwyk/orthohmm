# NumPy Binary Component Notice Mapping

This narrows one gap in the [wheel notice collection](INTEGRATED_DEPENDENCY_NOTICES_20260927.md).
Read the actual pinned NumPy 2.2.6 wheel and its embedded
`numpy-2.2.6.dist-info/LICENSE.txt`, not a current website or a different wheel.
The wheel digest was checked before/after reading; the notice length and
SHA256 matched the retained inventory. Nothing was installed or executed.

The [machine-readable mapping](numpy_binary_component_mapping_20260927.json)
records four provider declarations and three distinct shared-object payloads:

| Provider component | Notice line | Provider license label | File correspondence |
| --- | ---: | --- | --- |
| OpenBLAS | 76 | BSD-3-Clause | Literal pattern matches the bundled OpenBLAS shared object |
| LAPACK | 112 | BSD-3-Clause-Attribution | Same object; provider describes LAPACK as bundled in OpenBLAS |
| GCC runtime library | 167 | GPL-3.0-with-GCC-exception | No literal match; version-suffixed libgfortran filename candidate |
| libquadmath | 951 | LGPL-2.1-or-later | No literal match; version-suffixed filename candidate |

The provider's latter two patterns end in `.so`, whereas actual filenames end
in `.so.5.0.0` and `.so.0.0.0`. The receipt preserves literal failure and the
separate suffix-qualified candidates rather than silently widening a match.
All matched/candidate member bytes have length and SHA256 records. Four
declarations do not mean four separately bundled files.

The notice separately lists source-tree components and build tools it says
are not installed. Those statements must not be flattened into a claim that
all named components are independently identifiable wheel libraries. Likewise,
the inspected llvmlite third-party notice refers to LLVM source code; that
declaration alone is not a full compiled-component inventory.

This is provider-attribution bookkeeping, not legal advice or a clearance
determination. It does not establish exact component source versions, build
recipes, independent binary composition, static/transitive completeness,
source correspondence or satisfaction of redistribution obligations. All
license labels above are transcribed provider labels, not our interpretation.
No public runtime redistribution or publication-readiness claim follows.

To inspect the same evidence, use the wheel path and digest in the receipt,
read its recorded notice member with `zipfile.ZipFile`, verify the pinned
notice hash, and compare the four recorded `Files:` glob patterns against
the archive member names. The suffix candidates use the original pattern
plus `.*` and remain explicitly qualified. Hash each identified member's
uncompressed bytes; do not load its code.
