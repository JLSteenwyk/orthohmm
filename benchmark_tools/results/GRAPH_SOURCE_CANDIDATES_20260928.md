# Graph Dependency Source Candidates

The [acquisition receipt](graph_source_candidates_20260928.json) binds the
retained igraph and leidenalg wheels from the integrated workflow to their
exact PyPI release metadata and downloads the single source archive advertised
for each release. Both wheel SHA-256/size pairs match the published wheel
entries; both downloaded source archives match their published SHA-256/size.
Root PKG-INFO names and versions also match. No code was installed or executed.

| Package | Source archive bytes | Regular source files | Notice candidates |
| --- | ---: | ---: | ---: |
| igraph 1.0.0 | 5,077,105 | 2,909 | 15 |
| leidenalg 0.11.0 | 452,850 | 53 | 1 |

Source SHA-256:

- igraph: `2414d0be2e4d77ee5357807d100974b40f6082bb1bb71988ec46cfb6728651ee`
- leidenalg: `f454be96bbc8089ea2a90ca853d8d389ab646de964a03bd58417f8b29ff8ef5d`

Archives and raw PyPI responses are retained locally under
`benchmarks/work/graph_source_candidates_20260928`, not committed. The receipt
contains their identities and every regular source member's size and hash.
It is an acquisition manifest, not a source-to-binary build attestation.

## What This Adds

The igraph archive includes `vendor/source/igraph`, whose IGRAPH_VERSION file
contains `1.0.0`. Besides the wrapper LICENSE, it contains the C core COPYING
and notice candidates for bundled CSparse, f2c, GLPK, AMD/COLAMD, MiniSat,
Infomap, PCG and Qhull source directories. This provides source-level evidence
beyond the wheel's single top-level notice candidate. It does not establish
which optional code was linked into the retained wheel. The upstream
[installation documentation](https://python.igraph.org/en/main/install.html)
also distinguishes the Python wrapper from the vendored C core.

The leidenalg archive contains wrapper source and dependency-build scripts,
but not the two external library source trees those scripts download:

- `scripts/build_igraph.sh` selects igraph tag `1.0.0`.
- `scripts/build_libleidenalg.sh` selects libleidenalg tag `0.12.0`.

The retained wheel declares `libigraph-e2cb8a7d.so.4.0.0` and
`liblibleidenalg-32f15777.so.0.11.1`. Do not infer source versions from these
library names alone or treat the source archive's build scripts as proof of
the wheel's build recipe. For example, the tagged libleidenalg 0.12.0
[CMake target](https://github.com/vtraag/libleidenalg/blob/0.12.0/src/CMakeLists.txt)
sets SOVERSION 1 and VERSION from PROJECT_VERSION. The relationship between
the retained wheel and the external source revisions still needs evidence.
No library replacement, downgrade, or source-build attempt was made.

## Remaining Release Work

These are source candidates, not a complete corresponding-source bundle or
redistribution clearance. Exact external-library revisions, compiler/build
options, optional/static components and generated sources remain to be
established. In particular, the igraph wheel's bundled libxml2, liblzma and
libgomp are not accounted for merely by acquiring the igraph Python archive.
NumPy/Numba/LLVM components, external tools, Python and OS libraries still
require their own review. Filename-selected notices may miss inline notices
and are not a component-license mapping. No legal compatibility conclusion
is made here.

## Validation And Reproduction

72 focused source-acquisition, ELF-inventory and notice tests pass (0.53 s).
They cover exact wheel binding, source download mismatch, package identities,
ambiguous releases, unsafe names, archive links, duplicate members, redirects,
size limits and refusing overwrites. Inspection reads members without
extracting archive paths, importing modules or invoking build backends.
Failed downloads are left in their fresh output directory for inspection;
there is no silent resume or retry.

```bash
python -B -m benchmark_tools.acquire_wheel_sources \
  --inventory benchmark_tools/results/integrated_wheel_elf_20260927.json \
  --package igraph --package leidenalg \
  --output /fresh/graph-source-candidates \
  --receipt /fresh/graph-source-receipt.json
```

Scientific settings, scores, frozen runtimes and the timing panel are unchanged.
