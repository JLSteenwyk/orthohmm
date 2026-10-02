# Private Phylogeny Tools from Pinned Public Artifacts

MAFFT core and FastTree are now prepared from the independently acquired
public artifacts without historical tool installations or old fixture
inputs. All 34 MAFFT helpers match their historical identities. Only the
generated launcher's default helper-directory prefix differs. This removes
the old-installation prerequisite for private preparation on the recorded
host/compiler, not all runtime, numerical, scientific or release requirements.

## Controller and Tests

Checkpoint `7e6b12fd420f15cb379f31a43a70470d1208275f` adds
`prepare_publication_phylogeny_tools.py` and focused tests. Its embedded
34-helper catalog contains names/sizes/hashes, not source or binary payload.
That catalog derives from the independently pinned
[previous MAFFT build receipt](publication_mafft_build_20260926.json).
The old build driver and scientific admissions remain unchanged.

All **245 tests pass in 27.70s**, zero failures/errors/skips, across preparation,
acquisition, MAFFT/FastTree builds and source packaging. Retained JUnit is
37,762 bytes, SHA-256
`7e350d2df29f0f8f1bac697e19a4e662959f3d1b004681e9eb2c9428a4e0497d`.
Its suite time is 27.653 seconds; CLI wall time is diagnostic, not controlled
scientific efficiency. Earlier acquisition tests overlap, not additive.

Preflight binds an externally verified source index, seven exact public
artifact identities, Linux x86-64/AVX2 support and externally recorded GCC
digest. The AVX2 requirement applies to the retained upstream FastTree binary;
it does not establish a general requirement for OrthoHMM. Fresh canonical
output lies outside source/artifact roots. Shell metacharacters/spaces are
rejected because the upstream make recipes do not consistently quote PREFIX.
Acquired-file aliases/links, compiler/source changes and missing identities
fail before output creation.

Private HOME/TMP, fixed system PATH, empty make override variables and single
library-thread limits exclude inherited CC/LD_PRELOAD/MAFFT/Python overrides.
The core command uses `-j2`, explicit GCC, `CFLAGS=-O3` and the
fresh private installation prefix. It neither builds optional RNA extensions
nor changes source files or scientific settings. Two known absolute generated
convenience links are replaced with checked relative links; every installed
link must stay inside the prefix. Source/license/extension notices and all six
FastTree artifacts remain preserved privately.

Failures preserve logs, generated artifacts and `failed.json`; bounded
child-process groups stop on timeout with no retry. Tests inject each stage
failure, incorrect versions, changed helpers/modes/inventory/source/input,
unsafe links/paths, architecture/feature and compiler/source mismatches.

## Actual Copied-Source Preparation

Export/copy `native-build`: 1,844 payloads, consisting of 43 scientific,
1,800 workflow and one separate setup overlay. Total payload is 11,015,651
bytes; 1,823 Python files pass syntax checks. The 746,269-byte source index is:

`a9eac855a8677bc2bf35609f34c0a0bb3c11e441ac0d9cc7bc94db92abd6dd07`.

The four outer stages (export, copied verification before, private tool
preparation, copied verification after) return zero. The copied driver uses
isolated Python3.12 outside checkout; verification has an empty PATH. Actual
make/compiler/version child commands use the controlled system PATH. No
original tool prefix/source checkout/scientific Python package or old fixture
is an execution input.

Reuse all seven [acquired public artifacts](publication_phylogeny_acquisition_20261002.json).
No download is repeated, and every acquired input remains unchanged/non-executable.
The pinned source archive unpacks 173 original files; every original file
remains byte-identical after compilation. The recorded compiler is GCC13.3.0,
resolved executable `/usr/bin/x86_64-linux-gnu-gcc-13`, SHA-256:

`6117c52522997d2aaccb2b52b3c6bf42c0a6c5edb1d718431fed6b2fc5fec234`.

All five native stages succeed: compiler version, MAFFT core build, launcher
version, compiled helper version and FastTree help. Observed versions are
MAFFT `v7.525 (2024/Mar/13)`/helper `7.525` and FastTree2.2.0
double precision. These are version probes, not scientific inference.

All **34 helpers/13,744,200 bytes** match the historical sizes/hashes and
mode0755. The prepared tool inventory has **48 regular files and 18 in-prefix
symlinks**, with both known convenience links relative. FastTree's staged
copy becomes executable; the acquired original remains unchanged at mode0644.
Preserve the **147 warning occurrences** in the MAFFT build log, matching
the earlier build count; success is not a warning-free or memory-safety claim.

The generated launcher is 109,908 bytes, SHA-256
`558056349ca625e1dcfcfa8266be63663ea3523e6529c07dc9025be1c4a1c244`.
Independent byte-line comparison against the prior launcher confirms exactly
one changed line: line39's default helper-prefix assignment. Old launcher
bytes/receipt are read only after successful preparation; the controller
requires neither. The explicit `MAFFT_BINARIES` override selects the
new helper directory and must be adjusted if relocating the tools.

## Retained Evidence

The [machine execution receipt](publication_phylogeny_preparation_20261002.json)
is 133,716 bytes, SHA-256
`735177ca8c4f4705269628b947605ae0d4df884cf35d739585e0912f6e113036`.
Its copied helper completion is 101,744 bytes, SHA-256
`81f0fb18cf02e96c7976a56002580e5bfc36f866a9b367dccb72b8ed9fc85fe1`.

The [independent readback](publication_phylogeny_preparation_validation_20261002.json)
verifies 343 recorded file identities, 48 file-mode/18 link records, all stage
outcomes, 245-case JUnit, unchanged acquired inputs and the sole launcher-prefix
delta. It is 1,998 bytes, SHA-256
`428fb6bcaba2bb9e21f320fd1b87359fbd5bd8adb254ddd998325139ad7a4cc0`.
Repeated record references are integrity checks, not independent experiments.

Exact preceding [CI37048822738](https://github.com/JLSteenwyk/orthohmm/actions/runs/37048822738)
for pushed `5d049544` remains in progress at the retained
18:50:04UTC observation: docs, CPU wheel and native diagnostics are successful,
while five test jobs are in progress. Its
[status receipt](ci_phylogeny_preparation_prior_revision_20261002.json) does
not infer test counts, final CI success or CI status for this newer source.

## Remaining Boundary

[Reader commands](../PUBLICATION_SOURCE_COMPONENT.md#offline-private-phylogeny-tool-preparation)
use fresh paths and external anchors. Private preparation is not a hermetic
compiler/OS closure or universal reproducible build: different compiler
helper bytes fail the historical catalog. FastTree remains the exact retained
upstream binary, not the separate source-built baseline binary. No optional
RNA engine, full numerical alignment/tree check or full-dataset validation
is newly performed. Earlier native/scientific/collector admissions are reused
only within their existing bindings.

The fresh launcher is not silently substituted in the historical asset
manifest or admitted executor. Explicit source/wheel/base/tool assembly and
executor validation remain next handoff work. Compiler/OS and selected-file
security/rights closure, QfO uncertainty, independent-validation limitations,
matched resource evidence, final manuscript/release/archive/DOI and overall
publication readiness remain unfinished. No benchmark score/default/lock,
unrelated job/service, shared environment or DGX operation changes. Timing
remains deferred and the original publication goal stays active.
