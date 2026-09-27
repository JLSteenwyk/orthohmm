# Relocated Recovery Entrypoint

The new `export_recovery_entrypoint.py` exports three unchanged runtime harness
files and the dependency lock from committed revision `ebac5b2`. It refuses an
existing destination and verifies a relative file inventory with byte counts
and SHA256 values. Its manifest is change detection, not a signature. Wheels,
external tools and datasets are prerequisites, not included artifacts.

```bash
python -m benchmark_tools.export_recovery_entrypoint \
  --repo /absolute/repo --output /absolute/fresh/export
```

The [executed fixture receipt](relocated_entrypoint_fixture_20260927.json)
records a run from `/tmp/orthohmm-recovery-entrypoint-20260927`, with relocated
harness, input and output paths and isolated Python (`-I`). Its four input
files match the retained 16-gene fixture byte-for-byte. The existing validated
recovery interpreter and original MAFFT/FastTree paths were deliberately reused.

Native execution completed once with exit zero. All four independent readers
passed: three root groups, 36 ortholog pairs, one reconciled family, two bypassed
families, four duplications, three speciations and two species-tree marker
families. No species/gene-tree checkpoint was reused. The root partition is
byte-identical to the earlier full-pipeline fixture. All 13 focused exporter
and entrypoint tests pass, including mutation, unexpected inventory, symlink
and existing-destination rejection.

The local archive under `benchmarks/work/publication_relocated_entrypoint_20260927.tar.gz`
preserves the exported files, copied fixture inputs, outputs and execution/audit
receipts; its checksum is in the report. Large/raw artifacts were not added to
Git. This is not a published archive or redistribution clearance.

This establishes a small-fixture harness relocation, not a fresh environment
installation, complete original-path independence, full-data relocation or
cross-platform portability. No filesystem-access trace was collected. The
fixture has an empty satellite merge trace, and the independent readers still
ran from the original checkout. The diagnostic helper's other entrypoints
are not supported by this minimal export. A complete release still needs
environment/tool provisioning and independent execution on the target host.
