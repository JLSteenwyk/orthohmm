# Baseline FastTree Source Build Result

The [prespecified build protocol](FASTTREE_BASELINE_BUILD_PROTOCOL_20260926.md)
was committed in `9f9f6ee` before execution. Both private builds from pinned
FastTree 2.2.0 source completed with identical 413,192-byte binaries:

`e3f5f951d1b8e805eee5efd68b454f046a54be8a9eb91d3b8e2c6e5152aa9a30`

GCC 13.3.0 compiled double precision with explicit `-march=x86-64
-mtune=generic`, no native/AVX flags, OpenMP, single precision or fast-math.
Each build emitted three warnings, retained in its log; source was not patched
to suppress them. The binary's ELF note declares `x86-64-baseline`; observed
dynamic dependencies are libc, libm and the Linux loader. This is local
same-toolchain repeatability, not a hermetic or cross-platform build proof.

## Preserved Validator Failure

The [original driver receipt](publication_fasttree_baseline_build_failed_20260926.json)
remains failed. Its `-help` probe exited zero and printed the correct version,
but the driver incorrectly expected the word `Version` in that banner. The
pinned source explicitly prints the shorter banner for `-help`. The correction
and separate retained-binary verifier were committed in `a315c38` before the
first fixture execution. Neither build was repeated, no compiler flag changed,
and the failed report was not replaced.

## Installed-Pipeline Verification

The [separate verification receipt](publication_fasttree_baseline_verified_20260926.json)
confirms that the retained build runs with the rebuilt MAFFT and the installed
frozen-source OrthoHMM wheel on the unchanged 16-protein fixture. All five
prespecified outputs match the historical AVX2-backed fixture byte-for-byte:
orthogroups, root HOGs, native pairs, rooted species tree and reconciled gene
tree. All four independent readers also pass:

- [Tree/pair structure](publication_fasttree_baseline_structure_20260926.json).
- [Candidate residues, alignments and supermatrix](publication_fasttree_baseline_sequence_20260926.json).
- [Node events, final groups and complete native pairs](publication_fasttree_baseline_events_20260926.json).
- [Reconciliation selection and hierarchy](publication_fasttree_baseline_hierarchy_20260926.json).

All 123 focused build, acquisition and phylogeny-reader tests pass, including
banner regression tests and preservation of injected fixture failure. No
third-party source/binary was committed, and no old tool installation changed.
Job 22179 continues using its originally pinned upstream FastTree binary.

## Reproduction and Limits

For a new build directory with the corrected driver:

```bash
python -m benchmark_tools.build_publication_fasttree \
  --repo . --source benchmarks/work/publication_fasttree_source_20260926/upstream \
  --mafft benchmarks/work/publication_mafft_source_20260926/build/install/bin/mafft \
  --output benchmarks/work/publication_fasttree_baseline_new
```

The corrected driver's complete fresh-directory path has not been re-executed;
the actual retained-build path and its fixture/readbacks are recorded above.
The new binary is not admitted as a replacement publication benchmark runtime.
No full-dataset or older-CPU execution, speed advantage, statistical equivalence,
or complete dependency/redistribution clearance is claimed. The source header's
GPLv2-or-later declaration and separate GPLv3 LICENSE remain preserved locally;
compiler/OS and transitive source obligations still require artifact-level
review. This advances runtime packaging, not publication readiness by itself.
