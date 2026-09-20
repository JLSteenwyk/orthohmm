# Patched CPU Wheelhouse Notice Inventory

Inventoried the exact11local wheels named in the retained patched baseline
CPU installation report. Every wheel matched its recorded SHA-256, filename
distribution/version and embedded top-level metadata. This is the development
installation artifact set, not a new freeze of the scientific executors.
No dependency was installed, loaded, extracted, modified or redistributed.

| Distribution | Version | Provider metadata declaration | Notice candidates | Native members |
| --- | --- | --- | ---: | ---: |
| DendroPy | 5.1.0 | BSD | 3 | 0 |
| igraph | 1.0.0 | GNU General Public License (GPL) | 1 | 4 |
| leidenalg | 0.12.0 | GPL-3.0-or-later | 1 | 3 |
| llvmlite | 0.49.0 | BSD-2-Clause AND Apache-2.0 WITH LLVM-exception | 2 | 1 |
| numba | 0.67.0 | BSD | 2 | 14 |
| numpy | 2.2.6 | Long license prose; BSD classifier | 4 | 22 |
| orthohmm | 0.5.0 | No License/License-Expression header; embedded project license present | 1 | 3 |
| pip | 26.2.1 | MIT | 44 | 0 |
| python-igraph | 1.0.0 | GNU General Public License (GPL) | 1 | 0 |
| texttable | 1.7.0 | MIT | 1 | 0 |
| setuptools | 83.0.0 | MIT | 19 | 0 |

These are verbatim short metadata labels/expressions, not independent
interpretation or compatibility determinations. NumPy's55,523-character
License field is retained by length/hash rather than summarized as a single
grant. The OrthoHMM wheel's project license was already checked against the
repository MIT file in its separate wheel verification. Absence of a metadata
header is not absence of a license document.

## Evidence And Scope

[Machine-readable inventory](publication_dependency_notices_20260920.json)
SHA-256: `69e332c1b1f5d1a8b696ad3af975b20b4c05f7ff8177e2209dc26bdaa7948d4f`.
It records all11wheel identities, metadata identities,79notice candidates,
47native-library members, declared-license-file matches and source identity.
Every declared license file resolves to exactly one candidate in this set.
That does not establish that every relevant notice was declared or found.

Notice candidates include declared files, members under a `licenses`
directory, and filenames starting with conventional LICENSE/LICENCE,
COPYING, NOTICE, COPYRIGHT, PATENTS or AUTHORS tokens. Candidate names and
bytes are hashed, not treated as a complete legal interpretation. AUTHORS
may be attribution rather than a license. Native member detection is by
filename, not a linkage/static-component analysis.

The initial real inventory attempt rejected setuptools because it contains
13metadata files: one top-level distribution plus12vendored distributions.
The scanner now selects top-level wheel metadata and still inventories
notice candidates in vendored directories. Regression coverage preserves
rejection of multiple top-level metadata files. The failed attempt produced
no result file and modified no wheel. This inventory does not enumerate all
vendored component versions or bind every notice to every compiled component.

Provider declarations and bundled notices must remain distinct from the
project's own MIT license. In particular, the GPL-family declarations and
llvmlite third-party license document need artifact-specific review before
any combined binary/environment archive is represented as cleared for public
redistribution. No compatibility verdict or blanket clearance is supplied.
External binaries, OS libraries, containers and datasets remain outside scope.

## Reproduction

From the repository root, with the retained local wheelhouse and a fresh
output destination:

```bash
python -m benchmark_tools.inventory_dependency_notices \
  --install-report benchmarks/work/publication_baseline_patched_install/install_report.json \
  --report-sha 8546a5162a6250293c506b1d61b6b06faea5bcea11eb936ab0ae7e9d2c9c4de7 \
  --output /tmp/orthohmm-dependency-notices-new.json
```

The tool only reads report-pinned local `file:` URLs; it rejects remote
URLs, hash/name/version mismatches, duplicate members, unsafe declared
paths and duplicate distribution rows. It refuses an existing report
destination and rechecks input hashes after inventory. It does not validate
the current installed environment or substitute for a package security audit.

Seventeen inventory tests and four existing project-wheel tests pass
(21total). The retained inventory is a release-review input, not a public
release, archival deposit or proof that the broader publication goal is met.
