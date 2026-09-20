# Frozen Scientific Source Archive

Prepared a local source-only archive of scientific revision
`7f3a9e40dd7e79f842cc2c11fb8b548f9a802806`. Its own `version.py` remains
0.5.0. The revision-qualified archive name does not create a new package
version, Git tag, public release or deposition identifier.

Artifact:
`benchmarks/work/publication_frozen_source_20260920/orthohmm-0.5.0-scientific-7f3a9e4-mode0022.tar.gz`

- Size:126,207bytes, with551,096uncompressed source bytes in43regular files.
- SHA-256:`8e8e1ea4839fe0a78c3a6b7410951e65b8ee4311dd99b7a4b00d05275a51b8f4`.
- Selected Git paths: `LICENSE.md`, `README.md`, `requirements.txt`,
  `setup.py`, and the complete `orthohmm` package directory.
- Every regular member's bytes and permission bits match the frozen Git
  blobs. All32Python files compile syntactically without import/execution.
- A second export with the same explicit settings is byte-identical.

The archive includes project C/CUDA and experimental source under the
package, not compiled libraries, dependency wheels, benchmark datasets,
reference files, output predictions or competitor executables. It excludes
the rest of the repository, including tests, analysis workflows, manuscript
and figures. Historical README links can therefore refer to omitted files.
This is one source component of the future publication archive, not the
complete executable/reproducible study package.

## Verification And Reproduction

[Verification report](publication_frozen_source_archive_20260920.json)
records all43file paths, Git blob IDs, SHA-256 values, byte counts and modes,
the Python interpreter and verifier identity. The verifier does not extract
or execute source. It rejects missing/extra/duplicate members, source-byte or
mode changes, links, traversal paths and a pre-existing report destination.

From the repository root, choose new destinations:

```bash
git -c tar.umask=0022 archive --format=tar.gz \
  --prefix=orthohmm-0.5.0-scientific-7f3a9e4/ \
  --output=/tmp/orthohmm-frozen-source-new.tar.gz \
  7f3a9e40dd7e79f842cc2c11fb8b548f9a802806 \
  LICENSE.md README.md requirements.txt setup.py orthohmm

python -m benchmark_tools.verify_frozen_source_archive \
  --repo . --archive /tmp/orthohmm-frozen-source-new.tar.gz \
  --output /tmp/orthohmm-frozen-source-new.json
```

Git's initial default-mask exports were byte-identical to each other
(SHA-256:`a962c0faa3e90071ee6afbacb589cbd0dd4f9dfc87f6c544bc4c4b6d4fe28851`)
but failed strict mode verification because regular files were group-writable.
Both remain in the local artifact directory. New explicit-mask exports
preserve0644source modes. No frozen source or verifier condition was weakened
to admit the initial export. Repeat-byte equality is established on this Git/
compression implementation; semantic member verification remains the primary
cross-environment check.

Ten tests pass, including two real Git exports and corruption/expansion
rejection. Tests report22warnings from repeated syntax checks of two existing
invalid-escape string literals in the frozen parser/writer. The source is
preserved unchanged; syntax compilation is not a full interpreter-portability,
runtime-behavior or installed-build test.

## Remaining Release Gates

The historical setup still probes host-specific compiler flags, can skip
native kernels, and builds optional CUDA when available. Its dependency
declarations are not a complete pinned benchmark runtime. The later CPU-wheel
build controls are development changes, not silently substituted here. Do
not infer numerical equivalence, CPU-only installation or full phylogenetic
reproduction from this archive check.

The artifact remains local, reproducible from committed objects; it is not
newly uploaded or stored as a redundant binary in Git. Complete analysis/
runtime packaging, acquisition-only handling of restricted inputs, final
file-level rights/notice review, release metadata and external archival
steps remain incomplete. The publication objective is unchanged.
