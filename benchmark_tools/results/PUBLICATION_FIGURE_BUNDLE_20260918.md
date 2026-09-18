# Relocatable Figure Evidence Bundle

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
