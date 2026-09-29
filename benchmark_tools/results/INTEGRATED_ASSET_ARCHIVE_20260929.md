# Local Integrated Execution Asset Archive

Archived the selected execution assets of admitted job 22337 after validating
its pinned records with the controller's own path-preserving checker. This
local archive contains 98 regular files and 18 relative symlinks, plus 15
directories; regular payloads total 131,502,430 bytes. The uncompressed archive
is 131,604,480 bytes. [Restoration receipt and complete member inventory](integrated_execution_archive_restoration_20260929.json)
record archive identity, file hashes, modes and symlink targets.

The selected roots are `assets`, `readers`, `reader_wheels`,
`reader_requirements.txt` and `workflow.py`. Raw benchmark data, prior run
outputs, base Python, installer and OS libraries are excluded. This is a local
reproduction artifact, not redistribution-cleared material or a public release.

## Creation And Restoration

From the repository root, the archive was generated using GNU tar:

```sh
tar --format=gnu --sort=name --mtime=@0 --owner=0 --group=0 --numeric-owner \
  -cf benchmarks/work/integrated_execution_assets_20260929.tar \
  -C benchmarks/work/publication_integrated_full_ob_20260927 \
  assets readers reader_wheels reader_requirements.txt workflow.py
```

Before extraction, checked every archive member against the selected source
inventory, rejecting unexpected names/types, escaping links and differing
bytes or modes. All links are relative and resolve within the asset root.
Archive ownership/timestamps are normalized. Use the separately retained
archive checksum, verify the member inventory, and choose a fresh destination
before restoring; the creation command itself does not prevent overwrite.

The first restore's default umask changed FastTree from mode 0777 to 0775.
Its [failed mode check](integrated_archive_mode_failure_20260929.json) is retained.
A second, fresh restore used `tar --no-same-owner --same-permissions -xf ...`
and passes every member/hash/mode/link check. Historical mode 0777 is preserved
for this local exact-restoration experiment, not recommended as deployment
permissions. No executable was run from the first restored directory.

From `/tmp`, isolated stdlib-only Python loaded the restored controller and
validated the restored native assets, independent readers, both exact
environment locks and the newly acquired/rebound full OrthoBench data manifest.
All preflight checks pass. The successful root is
`/tmp/orthohmm-restored-execution-assets-20260929-v2`.

This closes local execution-asset archive restoration and preflight, not fresh
installation/inference using the restored archive, cross-host validation,
hermetic OS restoration or legal/security clearance. Historical absolute paths
remain in provenance manifests; no filesystem access trace was collected.
The next validation is execution from restored assets with independently
supplied base Python/installer and acquired inputs, without changing scientific
settings or replacing historical results.
