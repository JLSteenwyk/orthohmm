# Leiden Source Version And Build Evidence

The apparent difference between the build script's libleidenalg `0.12.0`
selection and the retained wheel's `liblibleidenalg-32f15777.so.0.11.1`
filename has a reproducible source-version explanation. It is **not** proof
of the wheel's exact corresponding source and does not explain the separate
high-CPM SIGSEGV.

## Source Mechanism

The upstream `0.12.0` lightweight tag currently identifies commit
`7ce12032cd61a38caebc52881ce1f1a1be4f0d28`, tree
`39b3bf92e831558a5d23ad2672370cf533aae31c`. At that exact commit:

```text
git describe COMMIT                     -> 0.11.1-13-g7ce1203
git describe --tags --exact-match COMMIT -> 0.12.0
```

Plain Git describe considers annotated tags by default, whereas `--tags`
also includes lightweight tags. [Git documentation](https://git-scm.com/docs/git-describe).
The retained older `0.11.1` tag is annotated. The
[version source](https://github.com/vtraag/libleidenalg/blob/7ce12032cd61a38caebc52881ce1f1a1be4f0d28/etc/cmake/version.cmake)
calls `git_describe(PACKAGE_VERSION)` without requesting lightweight tags.
There is no committed VERSION or NEXT_VERSION override. It then takes the
prefix before the first hyphen, which gives `0.11.1` in the observed history.
The [target configuration](https://github.com/vtraag/libleidenalg/blob/7ce12032cd61a38caebc52881ce1f1a1be4f0d28/src/CMakeLists.txt)
uses PROJECT_VERSION for the shared-library VERSION property.

This source interpretation is consistent with the wheel filename. CMake was
not executed and the library was not rebuilt, so it is not an observed build
product. In particular, the filename alone is not evidence that the wheel
contains the old `0.11.1` source rather than the selected `0.12.0` source.
The [offline audit](libleiden_source_version_20260928.json) records all 45
regular source blobs, their Git/SHA-256 identities, selected build/notice
texts, tag inventory and command results. Git object integrity checks passed.

## Upstream Publication Evidence

The [metadata receipt](graph_upstream_build_evidence_20260928.json) retains
identities for ten public responses and both exact retained wheels.
Leiden's successful [release workflow 18978271528](https://github.com/vtraag/leidenalg/actions/runs/18978271528)
ran from wrapper commit `c836ef9cd72374490e016c5841098d6eb911a700`, tag
`0.11.0`, on 31 October 2025. All nine listed jobs succeeded, including the
x86-64 wheel build and PyPI upload. Both dependency-build scripts fetched
from that immutable commit exactly match the already retained PyPI sdist
members. The workflow calls those scripts for Linux wheel builds.

No release artifacts were listed by the GitHub API at this observation time.
Both exact wheel queries to the [PyPI Integrity API](https://docs.pypi.org/api/integrity/)
returned HTTP 404 with an explicit no-provenance message. These bounded
observations do not establish that historical artifacts never existed or
that no evidence is available elsewhere. A successful workflow record and
matching scripts are not an artifact-digest chain linking the native source
to the retained wheel, and no signature verification claim is made.

## Preserved Source And Checks

A source-only clone, a 51,070-byte Git-export archive and a 2,856,754-byte
history bundle are retained locally, not installed or committed as raw data.
The [preservation receipt](libleiden_preserved_source_20260928.json) records
hashes and commands. The bundle includes the release commit and old annotated
tag history. A fresh offline repository fetched only from that bundle and
reproduced both version labels; see the
[independent replay receipt](libleiden_bundle_replay_20260928.json).

Git archive applies the committed export-subst attribute: its exported
vcpkg.json contains `0.11.1-13-g7ce1203`, unlike the literal placeholder in
the source blob. Preserve that distinction when comparing archives or
reconstructing a source-only build.

```bash
python -B -m benchmark_tools.audit_libleiden_version \
  --repository benchmarks/work/libleidenalg_source_0_12_0_20260928 \
  --output /tmp/libleiden-version-audit-fresh.json
```

63 focused version-audit, source-acquisition, ELF-inventory and notice-export
tests pass. Controls demonstrate the difference between annotated and
lightweight release tags, committed objects versus worktree edits, and
rejection of wrong origins/commits, version overrides, symlinks and changed
version derivation. This is source-provenance work, not inference validation.

Exact historical compiler/options and source-to-wheel correspondence remain
unresolved, as do transitive component and redistribution reviews. Current
tag metadata is not automatically historical tag metadata. No package,
library, scientific configuration, benchmark score or runtime pin changed.
