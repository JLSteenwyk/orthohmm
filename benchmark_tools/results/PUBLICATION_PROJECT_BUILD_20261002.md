# Frozen Source Build and Exact Project Wheel Reconstruction

The historical project wheel is now reproducible from the frozen scientific
source, its separately pinned setup overlay, public build dependencies and
the recorded local GCC toolchain. All 42 source-built payload members match
the admitted wheel; reconstructing its ZIP metadata restores the **exact**
144,444-byte archive SHA-256
`cfdfde5ed1be29e4080dd3571f5c0fc5fc5f45b57ebe5ee9b3c7559096c3b93d`.
The build/reconstruction do not need the original project wheel as an input.
The original is read only for the independent post-build payload comparison.

This removes the previously mandatory local project-wheel archive dependency
for this artifact on the recorded toolchain. It does not publish a wheel,
close compiler/native/OS/rights dependencies, establish universal byte
portability, rerun scientific benchmarks or complete publication readiness.

## Explicit Source Profile and Tests

Checkpoint `01249971` adds `native-build`, extending native-wheel support
with a separate `build/setup.py` from revision
`6fd6df19daba83ec6467b917988f99e27a95be14`. Its fixed SHA-256 is
`88120e9d722557337d466a5026c4b238d5a21c9cc213bf35c27d4b480b347193`.
Scientific `setup.py` and all package sources remain unmodified at
`7f3a9e40dd7e79f842cc2c11fb8b548f9a802806`. Previous source profiles and
their verification remain unchanged. No default, endpoint or score changes.

The offline builder verifies source/input identities before installation,
creates a private venv, checks exact build dependencies and their payloads,
requires all three ELF kernels, checks metadata/source members and actually
loads kernel symbols. There is no accepted compiler fallback or optional CUDA.
Baseline flags and fixed `SOURCE_DATE_EPOCH` are explicit.

Initial five-module tests pass 255. Corrected source at `fb519b61` fixes
installation targeting and guards the supplied base; one test-only misplaced
fixture assertion fails (259 pass/one fail). Retain that JUnit record.
Checkpoint `d52e7d93` moves only that assertion to its proper test.
All **260 tests pass in 28.56s**; final six-module build/reconstruction/source/
handoff/acquisition panel passes **271 in 29.88s**, zero failures/errors/skips.
Counts overlap; do not add them. These are CLI wall times, not controlled
scientific efficiency measurements.

## Failed Attempt and Private-Base Invalidation

The first copied execution (`01249971`) stops at dependency checking,
before compilation. Its pip bootstrap runs under Python3.10 `-S` without
an explicit installation prefix, so pip writes into the **disposable supplied
private base**, not the fresh build venv. Installation logs and a later
snapshot confirm exactly pip26.2.1 and added setuptools83.0.0 in that prefix;
the build venv has no pip. This is consistent with the
[documented pre-3.14 site-dependent venv prefix behavior](https://docs.python.org/3/library/sys.html#sys.prefix).

Keep the failed build, both failure receipts and the altered prefix. The
earlier [private-base execution](PUBLICATION_PRIVATE_BOOTSTRAP_20261002.md)
remains historical evidence of its then-successful installation; its base
prefix is now explicitly **invalid for pip-only admission**. Do not silently
reuse its earlier distribution inventory as current or clean it to hide this.
The 253,984-byte invalidated-state snapshot has SHA-256
`ff993fc3d15eb35d4dc7c8d853372bf10e8fbcc8617dedb37a5f8b496957455d`;
it contains 1,645 site files and the two distributions.

The correction explicitly [targets the build venv prefix](https://pip.pypa.io/en/stable/cli/pip_install/#cmdoption-prefix),
with `--ignore-installed`, and requires pip-only base metadata plus unchanged
site-file size/hashes before/after. Outputs inside the supplied base are
rejected. Provision one justified fresh private base using the existing
private Miniforge and retained public archives/pip; no artifact is downloaded
again. That base has all 19 expected package triples, CPython3.10.13, only
pip26.2.1, and 475 matching pip wheel payload files.

No shared installation, user environment, unrelated process, service or DGX
operation is involved. No expensive scientific/native fixture is restarted.

## Actual Corrected Build

Export/copy `d52e7d93` outside checkout and anchor its source index:
1,838 payloads (43 scientific, 1,794 workflow, one overlay), 10,942,341 bytes.
The 743,732-byte index SHA-256 is
`7e76abf0d388d64cfdcf9f95ede994a4cd68d8dcf15294ad2ab02ffde2689ac7`.
Isolated copied verification passes before/after, Git absent from PATH.

The build actually uses pinned public pip26.2.1/setuptools83.0.0 and system
GCC13.3.0. All nine controller stages pass, including library load and the
base-unchanged probe. All 883 base site files have unchanged identities,
including generated files; this is site scope, not every base/OS path.
Installed build payload checks match 475 pip and 440 setuptools files.
All 33 frozen scientific wheel members match their source Git blobs.

| Kernel | Bytes | SHA-256 |
| --- | ---: | --- |
| `hmm_viterbi.so` | 21,008 | `eeb4985e6f35689497a9c6187db56f8094c2cb9d92af18b77c027af0dfbea1d2` |
| `kmer_prefilter.so` | 20,632 | `4bef1865bccee6c4a5df6237880cf22b272cd567a1dc2d918815a28c86c7ce5c` |
| `pair_align.so` | 20,784 | `1a1036e00730d9548d49819cdd852dc6aca74ace53c27ca2c903720904bc1d02` |

Every kernel is byte-identical to its historical counterpart. The source-built
candidate wheel is 144,444 bytes, SHA-256
`2b4dfa20491ac3dcbad0d738d1632cf7a2add9330e7e6376f8a40441262d2d43`.
All **42/42** decoded archive members match, including distribution metadata
and RECORD; archive wrapper bytes differ. The [build receipt](publication_project_build_20261002.json)
is 72,684 bytes, SHA-256
`91261a6f28fe81445e12cbc58a8b7c7009b63d384e3cce3b5dae3e554f58df97`.
Seven actual preparation/build/verification stages return zero.

## Exact Wrapper Reconstruction

Checkpoint `1363fcc0` adds a metadata-only historical recipe: names,
member hashes/sizes, order, timestamps and attributes. It embeds no source,
compiled payload, original archive, credential or private workstation path.
Preflight rejects any changed, missing, extra, duplicated or oversized member
before creating output. Serialization restores only the wrapper; final
size/SHA must equal the historical artifact, or the partial/failed attempt is
retained without retry. Different compiler payloads are never replaced.

Actually export/copy the newer native-build source: 1,840 files, 10,960,685
bytes, with 1,819 Python syntax checks. Its 744,583-byte index has SHA-256
`dfd80f007445ea2a8d8a033b82b2477bf528e550a09d9be90fbf17c31f75a98e`.
The copied helper runs once on the existing candidate under CPython3.10.13/
zlib1.2.13. It produces the exact historical wheel and passes the existing
acquisition helper's unchanged filename/size/hash checks. No compilation,
installation, native kernel loading, network or inference is repeated here.
The old 1,838-file native-build component still verifies with its original
external digest under the current verifier.

All six reconstruction/export/verification stages return zero.
The [reconstruction receipt](publication_project_reconstruction_20261002.json)
is 16,720 bytes, SHA-256
`314f4b63e29edaac1addbaacad5e884b4f10e0f17e11d2cd5177332ee95d0a5a`.
The 2,177-byte helper completion receipt has SHA-256
`b0f0de1e59168537fec472fbfbcca74a88140017ecb23b9309d44ae3b4f3e901`.
[Executable reader commands](../PUBLICATION_SOURCE_COMPONENT.md#frozen-source-cpu-wheel-build)
use separately public build artifacts; there is no circular need to supply
the historical project wheel before building it.

## CI and Remaining Work

The preceding pushed `14054cd1` has completed
[CI run37041128334](https://github.com/JLSteenwyk/orthohmm/actions/runs/37041128334)
successfully: all eight observed jobs successful at 18:15:07 UTC. Its
[status receipt](ci_project_build_prior_revision_20261002.json) records IDs,
without inferred test counts or a claim about these newer changes.

Exact project-wheel recovery does not close GCC/OS/static-source attribution,
security, redistribution, external native-tool acquisition, complete executor
handoff or public compiled delivery. Remaining QfO uncertainty, matched
resource evidence, final manuscript/release/archive/DOI and overall publication
readiness remain unfinished. Timing stays deferred without renewed quiet-window
questions or contention polls. The original goal remains active.
