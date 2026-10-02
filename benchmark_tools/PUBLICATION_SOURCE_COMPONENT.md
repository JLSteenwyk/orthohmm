# Publication Source Components

This source-only component supplements, rather than replaces, the retained
scientific source archive and figure/runtime evidence packages. It is not the
complete executable study, an installed package, a public release or a DOI.

## Contents

- `scientific/`: `LICENSE.md`, `README.md`, `requirements.txt`, `setup.py` and
  the complete `orthohmm` package from scientific revision
  `7f3a9e40dd7e79f842cc2c11fb8b548f9a802806` (historical version 0.5.0).
- `workflow/`: all committed top-level `benchmark_tools/*.py`, all committed
  `tests/unit/*.py`, the project license and the two publication source and
  reproduction guides from the separately selected workflow revision.
- Optional `build/`: only the pinned setup overlay from
  `6fd6df19daba83ec6467b917988f99e27a95be14`, with `native-build`. The original
  scientific `setup.py` remains separately preserved.
- `SOURCE_INDEX.json`: exact Git revision/path/blob, SHA-256, byte count and
  mode for every exported file. Its digest must be retained outside the bundle
  before transfer; internal consistency alone is not authenticity.

The optional `orthobench-inputs` profile additionally exports exactly three
frozen JSON manifests used by the existing acquisition verifier and rebinder.
It checks their fixed scientific hashes, requires both helper scripts, and
records schema v2 and the explicit profile in the index. Default `source-only`
exports and historical v1 verification retain their original selection.

The `native-preparation` profile extends that selection with the exact
historical base reconstruction receipt and required offline installer helper
modules. New workflow revisions also include the public-artifact acquisition
and private Miniforge bootstrap helpers. It does not include package archives,
the pip wheel, a Conda bootstrap or native tools. The four result-document
exceptions are explicitly pinned.

The optional `native-wheels` profile adds three further fixed documents: the
admitted 12-wheel inventory and exact inference/reader hash locks. It requires
the wheel and base acquisition helpers in addition to native preparation.
Earlier profile selections and their historical verification remain unchanged.
This seven-document exception carries metadata/locks, not wheel payloads,
native predictions or raw sequence/reference content.

The `native-build` profile keeps those seven documents and adds the separately
revisioned setup-only overlay. Its SHA-256 is
`88120e9d722557337d466a5026c4b238d5a21c9cc213bf35c27d4b480b347193`.
New workflow revisions include the offline project builder and metadata-only
historical ZIP reconstruction helper. Earlier profile selections/verification
remain unchanged; neither an overlay nor a rebuilt artifact changes the frozen
scientific version, historical locks or admitted benchmark executor.

Dataset/reference files, test samples, native predictions, scoring outputs,
result receipts/plans, figures, manuscript assets, third-party source archives,
dependency wheels, binaries and runtime/OS images are excluded. This avoids
newly redistributing those materials; it does not resolve every source-file
attribution or clear an arbitrary future release. Both project license copies
are retained. Do not interpret their presence as an override of third-party
notices or terms.

The only result-file exception in `orthobench-inputs` is the three acquisition
support manifests. They carry filenames, hashes and historical provenance,
not FASTA sequence content or reference-group memberships. Their inclusion
does not authorize redistribution of separately acquired upstream data.
The documents are unmodified historical records and also contain earlier
plans, aggregate scores and provenance. Those fields are not a current
result-table replacement or commands to execute on their original paths.

The directories are separate intentionally. Do not merge their revision
identities, install workflow source over the frozen scientific package, or
treat the historical unpinned requirements as the validated publication
runtime. The benchmark-only ordering and dependency amendments remain in the
external execution plans; source packaging neither applies nor approves them.

## Build And Offline Verification

Run the builder from a clone with the two committed revisions, using a fresh
output directory and an explicit workflow commit:

```bash
python -B benchmark_tools/bundle_publication_source.py build \
  --repo . --revision WORKFLOW_COMMIT --output /absolute/fresh/source-component
```

After recording the emitted index digest externally, verification needs only
the relocated component and a compatible Python interpreter, not Git, the
original checkout, historical workstation paths or scientific dependencies:

```bash
python3 -I -B /relocated/source-component/workflow/benchmark_tools/bundle_publication_source.py \
  verify /relocated/source-component --manifest-sha256 RETAINED_INDEX_SHA256
```

The verifier checks the external index digest, byte identities, selected paths,
component mappings, file modes and exact inventory; it rejects symlinks,
escaping/duplicate/missing/extra paths and performs syntax compilation without
executing or importing exported Python files. Syntax success is not API,
dependency, numerical-equivalence or native reproduction validation.

## OrthoBench Input Preparation

Use the optional profile to provide the verifier's and rebinder's manifest
dependencies in the exported tree:

```bash
python -B benchmark_tools/bundle_publication_source.py build \
  --repo . --revision WORKFLOW_COMMIT --profile orthobench-inputs \
  --output /absolute/fresh/source-component
```

Verify its index as above before using its scripts. Then acquire the upstream
data separately at the frozen commit. These commands require Git/network
for acquisition and standard-library Python for verification/rebinding; use
fresh acquisition and output paths:

```bash
git clone --no-checkout https://github.com/davidemms/Open_Orthobench.git /absolute/fresh/orthobench
git -C /absolute/fresh/orthobench checkout --detach 872d6f30592ab5ff837224db16a514b3f2bb916a
python3 -I -S -B /relocated/source-component/workflow/benchmark_tools/verify_orthobench_acquisition.py \
  --checkout /absolute/fresh/orthobench \
  --results /relocated/source-component/workflow/benchmark_tools/results \
  --output /absolute/fresh/acquisition.json
python3 -I -S -B /relocated/source-component/workflow/benchmark_tools/rebind_orthobench_data.py \
  --manifest /relocated/source-component/workflow/benchmark_tools/results/integrated_orthobench_data_20260927.json \
  --sha256 fcae062525eec61a11de876ae798acffc0f2fe9614466c8ce339a6c214666061 \
  --acquisition /absolute/fresh/orthobench --output /absolute/fresh/rebound
```

The verifier compares 95 upstream files with their frozen Git blobs, including
93 benchmark inputs. The rebinder preserves the ordered 12 FASTA, 70 reference
and 11 low-certainty identities but replaces their historical paths with local
paths. Supply `rebound/data.json` and its digest from `rebound/rebind.json` to
the separately provisioned integrated workflow. References are scoring inputs,
not inference inputs. No upstream scorer or native inference runs here.
The native wheel sets, aligner/tree builder, base Python and installer still
require separate preparation; this profile is not the complete executable study.

## Private Historical Conda Bootstrap

New workflow revisions provide `prepare_publication_bootstrap.py`, a stdlib
helper for Linux x86-64. Verify the exported source index before using it.
The fixed historical official Miniforge release is 25.3.1-0, not the latest
release or a current-security recommendation. The installer is 93,870,801
bytes, SHA-256
`376b160ed8130820db0ab0f3826ac1fc85923647f75c1b8231166e3d559ab768`;
the 104-byte provider checksum sidecar is separately pinned. Acquisition
checks both identities before any installation:

```bash
python3 -I -S -B /relocated/source-component/workflow/benchmark_tools/prepare_publication_bootstrap.py acquire \
  --output /absolute/fresh/bootstrap-artifacts --acknowledge-historical-runtime
python3 -I -S -B /relocated/source-component/workflow/benchmark_tools/prepare_publication_bootstrap.py install \
  --installer /absolute/fresh/bootstrap-artifacts/Miniforge3-25.3.1-0-Linux-x86_64.sh \
  --output /absolute/fresh/private-bootstrap --acknowledge-historical-runtime
```

Each output must be fresh, canonical and outside the immutable source component.
The acquisition permits only fixed GitHub HTTPS URLs and release-asset redirects,
bounds size, checks hashes and excludes temporary signed redirect queries from
receipts. Socket timeouts are not a whole acquisition deadline. Failures retain
partial files/records without retry or overwriting an existing attempt.

Installation verifies the exact supplied artifact first, runs batch mode in a
fresh prefix with private HOME/cache/config and does not request `conda init`.
The private prefix's `bin` is first on PATH so an `env python` Conda launcher
resolves its own interpreter. It then checks Conda 25.3.1, internal entrypoint
paths and package metadata, retaining stage logs and `complete.json`. No shared
environment, shell initialization, service or scientific configuration change
is requested. This is downloaded installer execution after identity checks,
not a syscall sandbox or signed-binary/security/rights clearance.

Use `private-bootstrap/prefix/bin/conda` with its recorded entrypoint SHA-256
from `private-bootstrap/complete.json` for the offline base installer below.
That digest binds the installed entrypoint, not its interpreter or complete
bootstrap payload closure. Do not hardcode another installation's digest.
The [executed private-bootstrap check](results/PUBLICATION_PRIVATE_BOOTSTRAP_20261002.md)
records actual public acquisition, isolated installation and compatibility
with the unchanged frozen base. Full bootstrap/base/OS/rights closure and
public distribution remain separate requirements.

## Offline Historical Base Installation

Use `--profile native-preparation` with the builder to include the base
reconstruction document and controller. After verifying the emitted index,
the newer acquisition helper can download the 19 exact archives and pinned pip
wheel into a fresh private directory outside the immutable component:

```bash
python3 -I -S -B /relocated/source-component/workflow/benchmark_tools/acquire_publication_base.py \
  --receipt /relocated/source-component/workflow/benchmark_tools/results/reconstructed_base_fixture_20260927.json \
  --output /absolute/fresh/base-artifacts --acknowledge-historical-runtime
```

The helper requires the exact retained reconstruction receipt, uses only
credential-free allowlisted HTTPS endpoints/redirects and bounds each response.
All 19 archives must match their byte counts, SHA-256 and MD5. The release-
specific PyPI metadata must select the one pip wheel with the frozen size/hash;
no latest-version resolution or replacements. Verified downloads, metadata and
source/input pins appear in `complete.json`. Failures preserve records/partial
files and `failed.json` without success, retry or overwriting an old directory.
There is no package installation, shared-environment change or native inference.
Its timeout is a socket-operation bound, not a whole acquisition deadline.
HTTPS and recorded hashes do not establish security or redistribution clearance.

Pass `base-artifacts/archives` and
`base-artifacts/bootstrap_wheels/pip-26.2.1-py3-none-any.whl` to the offline
installer below. Alternatively, separately
supply a directory containing the 19 exact archives listed under
`acquisition.packages` in the exported
`results/reconstructed_base_fixture_20260927.json`. The document retains
their HTTPS provider URLs, SHA-256/MD5 and byte counts. Acquire those materials
separately; the controller neither downloads nor solves for replacements.

Supply the exact `pip-26.2.1-py3-none-any.whl` bootstrap artifact (1,816,632
bytes, SHA-256
`71138adf1f4ca900cdb7d289c21b7494329f2332b6d85f0e1c42108c0384ed3e`)
and a separately trusted, compatible Conda executable with its externally
recorded entrypoint digest. The preceding private-bootstrap helper is one
executed route; the offline base installer itself does not acquire Conda or
authenticate its entire transitive stack. Set `RETAINED_CONDA_ENTRYPOINT_SHA256`
to the `conda.sha256` value from that installation's verified `complete.json`
and substitute `/absolute/fresh/private-bootstrap/prefix/bin/conda` below.

```bash
python3 -I -S -B /relocated/source-component/workflow/benchmark_tools/install_publication_base.py install \
  --receipt /relocated/source-component/workflow/benchmark_tools/results/reconstructed_base_fixture_20260927.json \
  --cache /absolute/acquired-base-archives \
  --conda /absolute/trusted-conda/bin/conda --conda-sha256 RETAINED_CONDA_ENTRYPOINT_SHA256 \
  --pip-wheel /absolute/acquired-wheels/pip-26.2.1-py3-none-any.whl \
  --output /absolute/fresh/base-preparation --acknowledge-historical-runtime
```

Use an output directory outside the immutable source component. Complete
input preflight precedes creation; execution has private HOME/cache/config
and the supplied Conda entrypoint's directory first on PATH,
offline copied explicit packages, a hash-required pip overlay and no retry.
The four stages install Conda packages, bootstrap pip, run `pip check` and
snapshot the installed Python. Success additionally requires all 19 expected
name/version/build triples, Python 3.10.13 on x86-64, exactly pip 26.2.1 in
site-packages, byte-checked pip payloads, unchanged supplied/staged artifacts
and a retained `complete.json`. Failures preserve their partial directory and
logs, do not emit success, and must not be silently resumed or overwritten.

The resulting `base-preparation/python-runtime/bin/python` can serve the
separately provisioned integrated workflow's `--base-python` and
`--installer-python`. This is the historical reconstruction policy already
used for the admitted full OrthoBench reproduction, not a security-cleared
installation recommendation. Conda metadata equality does not verify every
prefix-transformed payload, and shared OS libraries remain outside scope.
This step does not install the scientific environments, infer orthology,
admit timing, or complete the full study release.

## Frozen-Source CPU Wheel Build

Export `--profile native-build` and verify its external source-index digest.
Use a separately provisioned, unchanged historical CPython 3.10.13 Linux x86-64
base with exactly pip 26.2.1, recording its interpreter digest externally.
The preceding private-bootstrap/base commands provide one route. Outputs must
be fresh and outside both the immutable source component and the supplied base.
System GCC is required in `/usr/bin:/bin`; that PATH must not expose `nvcc`.
Its executable/version is recorded, not supplied or fully dependency-closed.

The build needs only two supplied public artifacts, not an existing OrthoHMM
wheel. Pip is available through the base-acquisition helper. Obtain the exact
[setuptools 83.0.0 wheel](https://files.pythonhosted.org/packages/5d/40/e1e72872c6354b306daef1703549e8e83b4d43cfea356311bf722a043752/setuptools-83.0.0-py3-none-any.whl)
separately: 1,008,090 bytes, SHA-256
`29b23c360f22f414dc7336bb39178cc7bcbf6021ed2733cde173f09dba19abb3`.
The controller requires both exact filenames/sizes/hashes before installation:

```bash
python3 -I -S -B /relocated/source-component/workflow/benchmark_tools/build_publication_project_wheel.py \
  --component /relocated/source-component --manifest-sha256 RETAINED_INDEX_SHA256 \
  --base-python /absolute/base-preparation/python-runtime/bin/python \
  --base-python-sha256 RETAINED_BASE_INTERPRETER_SHA256 \
  --pip-wheel /absolute/base-artifacts/bootstrap_wheels/pip-26.2.1-py3-none-any.whl \
  --setuptools-wheel /absolute/acquired/setuptools-83.0.0-py3-none-any.whl \
  --output /absolute/fresh/project-build --acknowledge-historical-runtime
```

The controller copies scientific sources, replaces only staged `setup.py`, and
creates a private build venv without inherited site packages. Hash-required
installation explicitly targets that venv's prefix and ignores existing
installations; it does not rely on Python3.10 `-S` preserving a venv prefix.
It checks the two installed distributions and their wheel payloads, then builds
offline without dependency resolution/build isolation/cache, with baseline CPU
flags and fixed `SOURCE_DATE_EPOCH`.

Success requires exactly one platform wheel, all 33 frozen scientific members,
all three ELF kernels, expected distribution/tag metadata and actual kernel
loading/exported symbols, including a non-AVX2 baseline witness. Base site-package
metadata and every site file's size/hash must remain identical before/after;
source/staged inputs and the compiler alias are also checked. All nine stages
and their logs are retained. Failures preserve artifacts and `failed.json`
without retry or success. This is not a syscall sandbox or full base/OS audit.

A new wheel is not silently accepted by historical locks. The metadata-only
reconstruction helper can restore the original ZIP wrapper only when all 42
payload members already match their historical identities:

```bash
/absolute/fresh/project-build/venv/bin/python -I -S -B \
  /relocated/source-component/workflow/benchmark_tools/reconstruct_publication_project_wheel.py \
  --candidate /absolute/fresh/project-build/wheels/orthohmm-0.5.0-cp310-cp310-linux_x86_64.whl \
  --output /absolute/fresh/reconstructed-project-wheel
```

It embeds names, hashes, sizes, order, timestamps and attributes, not payloads.
It neither needs the old wheel nor recompiles, installs or loads project kernels.
Success requires the exact 144,444-byte historical archive SHA-256
`cfdfde5ed1be29e4080dd3571f5c0fc5fc5f45b57ebe5ee9b3c7559096c3b93d`.
Compression/compiler differences fail closed with retained attempts; another
compiler's different binaries are not substituted. Use the prescribed
Python/zlib/toolchain rather than assuming universal byte reproduction.
The reconstructed artifact can serve the historical acquisition controller's
supplied project-wheel input. This is an artifact-recovery route, not a public
PyPI release, complete runtime/rights closure or independent new accuracy result.
See [actual build, failures and reconstruction](results/PUBLICATION_PROJECT_BUILD_20261002.md).

## Historical Inference and Reader Wheel Acquisition

Export `--profile native-wheels` instead of `native-preparation` and verify
its externally anchored index before using the following acquisition workflow.
Supply the exact unpublished OrthoHMM setup-overlay wheel from the preceding
source-build/reconstruction route or the retained private execution archive,
144,444 bytes, SHA-256
`cfdfde5ed1be29e4080dd3571f5c0fc5fc5f45b57ebe5ee9b3c7559096c3b93d`.
The matching pip wheel can come from the preceding public base acquisition.
Do not treat OrthoHMM 0.5.0 as an available PyPI release or silently replace it
with an arbitrary local build. Without that exact supplied wheel, this workflow
does not produce the complete historical inference environment.

```bash
python3 -I -S -B /relocated/source-component/workflow/benchmark_tools/acquire_publication_wheels.py \
  --inventory /relocated/source-component/workflow/benchmark_tools/results/integrated_wheel_elf_20260927.json \
  --inference-lock /relocated/source-component/workflow/benchmark_tools/results/publication_recovery_requirements_20260926.txt \
  --reader-lock /relocated/source-component/workflow/benchmark_tools/results/publication_reader_requirements_20260927_v2.txt \
  --orthohmm-wheel /absolute/supplied/orthohmm-0.5.0-cp310-cp310-linux_x86_64.whl \
  --pip-wheel /absolute/base-artifacts/bootstrap_wheels/pip-26.2.1-py3-none-any.whl \
  --output /absolute/fresh/native-wheels --acknowledge-historical-runtime
```

The controller reads the exact retained inventory and locks, not the historical
wheel paths embedded in its metadata. It verifies the two supplied artifacts,
downloads ten exact public third-party wheels using release-specific PyPI
metadata and preserves names/platforms/sizes/hashes without solving or choosing
latest versions. Copies prepare 11 `inference_wheels` and five `reader_wheels`
with the unchanged locks at `inference_requirements.txt` and
`reader_requirements.txt`; the union has 12 artifacts. The prior safe HTTPS
and response bounds apply, with partial/failure retention and no retry.

The inference set contains historical Leiden 0.11.0, not the different 0.12
overlay lock. The reader set contains the admitted Biopython 1.87 amendment,
not the retained 1.86 history. `complete.json` binds all supplied/downloaded
wheel bytes, provider metadata, copies and exact locks. This prepares artifacts,
not an installed/compatible/security-cleared runtime, scientific source audit,
declared dependency resolution or inference admission. Installation, tool and
reader-source assets, bootstrap/OS and rights obligations remain separate.

## Executable Study Boundary

The workflow reproduction guide is included for context, but links to omitted
result evidence and figures remain repository-relative and are not rewritten.
Native execution still requires the pinned dependencies, acquisition-only raw
inputs, exact command/data manifests and external tool/runtime assets described
there. Some historical scripts contain DGX commands and absolute paths: their
inclusion preserves history, not authorization to run them. The approved timing
host remains the local Threadripper; no DGX operation is required or authorized.

Unit-test source is included but test execution can require excluded fixtures,
dependency installations, native tools or retained receipts. This component
does not assert that every test can run in isolation. Scientific conclusions,
failed analyses, controlled timing admission, redistribution review and final
publication completion remain governed by the retained evidence and goal.
