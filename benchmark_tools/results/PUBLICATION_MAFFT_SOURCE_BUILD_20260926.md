# MAFFT Source Build and Installed Pipeline Check

Reacquired the official MAFFT 7.525 with-extensions source archive over HTTPS,
SHA-256 `2876f4adc1a2de4ed206bc40896763bf208bf1a02bda52f8bfdd91cf52d73e4a`.
Of 173 regular archive files, 172 match the retained installation's source
bytes. The sole difference is the first line of `core/Makefile`: the retained
copy replaces `/usr/local` with its installation prefix. No source was patched
for the new build; the new prefix is a make command-line assignment.

## Observed Results

- A private core-only build with GCC 13.3.0, `CFLAGS=-O3` and `make -j2 install`
  completed under a new prefix. No extension build, root access, historical
  installation change or DGX access occurred.
- All 34 installed helper files, including the two Perl helper scripts, are
  byte-identical to their historical counterparts. The launcher's embedded
  prefix differs by design. This demonstrates same-host reconstruction of
  these files, not an isolated compiler/toolchain or cross-platform build.
- The installed frozen-source OrthoHMM pipeline completes the same 16-protein,
  four-taxon fixture with fresh inferred phylogeny and 36 cross-species pairs.
  The partition, root-HOG table, pair table, inferred rooted species tree and
  reconciled gene tree are byte-identical to the original installed fixture.
- The readback audit verifies 719 record entries. Nineteen related tests pass,
  including archive path/link/duplicate/hash rejection and content comparison.
- The compiler log contains 147 `warning:` occurrences. These remain in the
  recorded log; success is not a warning-free-build or memory-safety claim.

Evidence:

| Record | SHA-256 |
| --- | --- |
| [Build](publication_mafft_build_20260926.json) | `4c4e92a29c1dc4b27e4b47f4010e9ea9648c17beb02b5a00995a59535394ac76` |
| [Installed fixture](publication_mafft_build_phylogeny_20260926.json) | `47c2229dd4bac0a7c0651a297cbf9b235f647377c26393e7b34316c612d16e5d` |
| [Independent readback](publication_mafft_build_readback_20260926.json) | `efdee0c8856a190b23cdca37222de21fddcfd756d68cf4c8f0ffc27ba97f318c` |

The source and build artifacts remain local under
`benchmarks/work/publication_mafft_source_20260926/`. Only scripts, tests and
identity reports are newly committed, not third-party source or binaries.

## Reproduction

Use fresh download/build/readback paths; existing build destinations are
refused. The build script verifies the archive hash before unpacking, rejects
links, special files, unsafe paths and duplicates, and imposes a size limit.
The build and inference driver have process-group timeouts with no retries.
Failures retain logs and a report. The installed-pipeline verifier additionally
checks the frozen installed-source evidence before inference.

```sh
mkdir -p benchmarks/work/publication_mafft_source_20260926
curl --fail --location --proto '=https' --proto-redir '=https' --max-time 90 \
  --output benchmarks/work/publication_mafft_source_20260926/mafft-7.525-with-extensions-src.tgz \
  https://mafft.cbrc.jp/alignment/software/mafft-7.525-with-extensions-src.tgz
/usr/bin/python3 -S -m benchmark_tools.build_publication_mafft --repo . \
  --archive benchmarks/work/publication_mafft_source_20260926/mafft-7.525-with-extensions-src.tgz \
  --installed /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/mafft-7.525-with-extensions \
  --output benchmarks/work/publication_mafft_source_20260926/build
/usr/bin/python3 -S -m benchmark_tools.audit_publication_mafft_build --repo . \
  --build benchmarks/work/publication_mafft_source_20260926/build/report.json \
  --output benchmarks/work/publication_mafft_source_20260926/readback.json
python -m pytest -q tests/unit/test_build_publication_mafft.py \
  tests/unit/test_audit_publication_mafft_build.py \
  tests/unit/test_verify_frozen_phylogeny_install.py
```

The build runs with `/usr/bin:/bin` PATH and explicit GCC, prefix and CFLAGS.
This is not a hermetic environment: system utilities/libraries and other
environment variables remain available. The phylogeny driver explicitly sets
`MAFFT_BINARIES` to the new helper directory and retains the existing FastTree
executable and installed Python package. Do not transfer full-dataset scores
to a new environment solely from this synthetic check.

## Extension and Notice Scope

The official [installation instructions](https://mafft.ddbj.nig.ac.jp/alignment/software/source.html)
separate the core build from optional RNA structural extensions. All historical
helper files were reproduced by the core target alone. Wrapper programs such
as `mccaskillwrap` are core build products; their presence does not establish
that external RNA engines were installed or tested. RNA extension modes are
not needed by the protein fixture and were not evaluated.

The source archive still contains extension source and its notices, including
separate conditions for Vienna RNA, MXSCARNA and ProbConsRNA. A core-only
build does not remove those obligations from redistribution of the whole
source archive. Prefer acquisition instructions while final selected-artifact
notice and dependency review remains open. This is not legal clearance for
all bundled source, helper scripts, OS libraries or optional engines.

No scientific defaults or retained benchmark scores changed. Controlled
scaling, remaining uncertainty/provenance analyses, final runtime packaging
and release review remain incomplete.
