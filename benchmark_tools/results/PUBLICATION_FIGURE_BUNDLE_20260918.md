# Relocatable Figure Evidence Bundle

## Eighteen-Panel Refresh

The refreshed inventory adds corrected-QfO search coverage and SwissTrees
sequence-control uncertainty without changing the original 16-panel archive.
All 18 panels passed the direct-byte audit: 61 output records and 101 distinct
directly referenced files. The committed export contains 122 files totalling
23,074,148 bytes, excluding its 73,200-byte manifest.

Source revision: `f2749d78d32b938007c7e7af6da3e2f90038e2fd`.
The exporter reads committed Git blobs; historical manifests remain unchanged.
The detached method helper retains the same frozen revision described below.

- [Inventory audit](publication_figure_integrity_20260918_v3.json).
- [Bundle manifest](publication_figure_bundle_20260918_v2.json), SHA-256
  `fb60aad87ea93d38d0f9793ed2e36fb1b739b13b1967b716a59d97ce7d8cf1d9`.
- Local archive `benchmarks/work/publication_figure_evidence_20260918_f2749d7.tar.gz`:
  5,884,001 bytes, SHA-256
  `3ea3e6cc6a83260a10e9584a28513e1fc1aa2cbb98932399d80223c6544534bd`.
- [Relocation verification](publication_figure_bundle_relocation_20260918_v2.json):
  the archive was extracted to `/tmp/orthohmm-figure-refresh.mPCYCB` and the
  bundled verifier ran there under `/usr/bin/python3 -I`. All 122 files and
  dependency mappings passed without project imports or historical-path reads.

Twenty-eight focused inventory/bundler tests pass, including explicit CLI
selection of a newly committed audit. The old default audit is preserved,
so earlier reproduction commands retain their meaning. Reproduce this refresh
using fresh destination paths:

```bash
python benchmark_tools/bundle_publication_figures.py build \
  --repo . --revision f2749d78d32b938007c7e7af6da3e2f90038e2fd \
  --audit benchmark_tools/results/publication_figure_integrity_20260918_v3.json \
  --output /tmp/orthohmm-figure-evidence-refresh
python3 -I /tmp/orthohmm-figure-evidence-refresh/benchmark_tools/bundle_publication_figures.py \
  verify /tmp/orthohmm-figure-evidence-refresh
```

This verifies retained evidence, not regenerated statistics, plots, native
inference or scoring. Transitive dependencies, raw datasets, portable native
execution and redistribution clearance remain incomplete. This archive is
local, not externally deposited or assigned a DOI. The newer statistical
workflow still needs its own isolated reproduction check. Earlier limits
and evidence below remain applicable; no timing or accuracy claim is upgraded.

## Scope

All 16 explicitly retained figure panels now have a relocatable direct-evidence
bundle: 55 output records, 91 distinct directly referenced files, 16 original
figure manifests, the committed audit, repository license and standalone
exporter/verifier. In total, 110 files contain 22,068,376 bytes, excluding the
66,233-byte bundle manifest. Original figure manifests are preserved unchanged.
Their historical absolute paths are mapped explicitly to relative bundle paths.

This is an evidence-preservation component of the publication package, not the
complete publication archive. The verifier does not regenerate statistics,
figures, searches, trees, predictions or benchmark scores. Transitive imports,
compiled environments and raw data remain outside this bundle. Diagnostic
simulation results and descriptive DGX timing observations retain their original
limitations; packaging does not admit them as controlled comparisons.
Third-party redistribution clearance remains unfinished. The archive is local;
it has not been deposited externally or assigned a DOI.

## Provenance And Verification

Source revision: `fcb41ed094d3bc24329175c59b14c92c722e1d25`.
Every exported file comes from a committed Git blob, not the mutable checkout.
The one detached method helper comes from frozen revision
`7f3a9e40dd7e79f842cc2c11fb8b548f9a802806`, with its original SHA-256 verified.
All other direct dependencies must match the byte counts and SHA-256 values
in the committed figure audit and original manifests.

The machine-readable bundle manifest is retained as
`publication_figure_bundle_20260918.json`. SHA-256:
`e168560edaa842a9a157cf6259e291d7c903d7cde1930e4f5c69bccbc3fa1289`.
It records each file's Git revision/path and every historical-path relocation.

Local archive:
`benchmarks/work/publication_figure_evidence_20260918_fcb41ed.tar.gz`.
Size: 5,452,374 bytes. SHA-256:
`aa8518c9ce0e46200ef891336bd3b91d78b6f09dbee157e366ad43c51bbc7943`.
The archive is not committed as a duplicate binary dataset.

The actual archive was extracted to `/tmp/orthohmm-figure-relocation.u9etV9`
and verified using `/usr/bin/python3 -I`, from that directory, with the bundled
standard-library-only verifier. All 110 files and all historical dependency
mappings passed. No Git command, project import or historical source path is
needed for the verification step. This is a same-host relocation check, not
cross-platform execution validation.

35 focused tests pass, including a synthetic bundle relocated after its source
repository was moved away, dirty-checkout independence, changed committed
evidence, corrupt/missing/extra files, symlinks, path traversal, duplicate
entries, incomplete mappings and existing-output refusal.

The manifest and archive hashes above are external integrity anchors. The
internal verifier alone is not an authenticity check against coordinated edits
to both files and their manifest. Compare the retained manifest hash before
trusting a transferred bundle, and the archive hash before extraction.

## Reproduce

From a clone containing the source revision and frozen helper revision, use
new destination paths:

```bash
python benchmark_tools/bundle_publication_figures.py build --repo . --revision fcb41ed094d3bc24329175c59b14c92c722e1d25 --output /tmp/orthohmm-figure-evidence
tar --sort=name --mtime=@0 --owner=0 --group=0 --numeric-owner --format=ustar -czf /tmp/orthohmm-figure-evidence.tar.gz -C /tmp/orthohmm-figure-evidence .
python3 -I /tmp/orthohmm-figure-evidence/benchmark_tools/bundle_publication_figures.py verify /tmp/orthohmm-figure-evidence
```

The exporter refuses existing output directories. Tar creation is shown for
a new archive path and should not overwrite an existing archive. Packaging
metadata is normalized; compressed archive byte identity across different
tar/gzip versions has not been tested. File identities and the bundle manifest
are the portable content checks. Regenerating scientific results remains a
separate requirement, partially addressed by the existing relocated SwissTrees
workflow and not established for all figures by this export.
