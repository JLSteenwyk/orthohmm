# Canonical QfO Python Source Archive

The completed assessment plan references 762 distinct repository-contained
Python files (1,502 entries before deduplication). One pinned source,
`benchmark_tools/assemble_ob_provenance_register.py`, subsequently changed
when historical OrthoBench metadata was incorporated. Rechecking the original
plan against today's working tree therefore correctly rejects that path;
this does not retroactively change the admitted execution.

A new source-only local archive preserves all 762 exact pinned files,
5,429,985 uncompressed bytes. The changed module was recovered from Git revision
`97a277afe4b2e28a9d1a10bbc8dd905415daccb8` and checked against the original
plan's size and SHA256, not accepted merely because a commit contained it.
All other members matched retained bytes. No original path, plan, execution
receipt or scientific output was rewritten.

Artifact:
`benchmarks/work/qfo_canonical_python_sources_20260927/qfo-canonical-assessment-python-sources.tar.gz`

- Size: 1,268,109 bytes.
- SHA256: `4951aa49982441e51da980a6e5c9db019137cfc34ba14748b33901a6f1ee0c51`.
- All archive members were independently read back against original pins.
- A second export is byte-identical.
- Nine tests cover deduplication, conflicting pins, missing/extra/changed
  members, duplicate members, symbolic links and traversal-member rejection.

The [compact receipt](qfo_canonical_python_source_archive_20260927.json)
binds the archive, original plan, full local manifests, recovery revision and
exporting source. This is preservation evidence, not a new scorer execution.
The current admission checker still requires original source paths; it has
not been weakened to ignore edits or silently substitute archived files.

```bash
python -m benchmark_tools.archive_qfo_assessment_sources --repo . \
  --output /absolute/new/qfo-python-source-archive
```

This component contains Python source bytes only, with normalized 0644 modes.
It excludes raw data, references, containers, executables, wheels, non-Python
configuration and other runtime files. It is not a complete executable archive,
historical permission attestation, rights-cleared redistribution or deposited
release. Original assessment accuracy/provenance admission remains separate.
