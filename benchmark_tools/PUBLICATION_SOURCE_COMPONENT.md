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
modules. It does not include package archives, the pip wheel, a Conda bootstrap
or native tools. The four result-document exceptions are explicitly pinned.

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

## Offline Historical Base Installation

Use `--profile native-preparation` with the builder to include the base
reconstruction document and controller. After verifying the emitted index,
supply a directory containing the 19 exact archives listed under
`acquisition.packages` in the exported
`results/reconstructed_base_fixture_20260927.json`. The document retains
their HTTPS provider URLs, SHA-256/MD5 and byte counts. Acquire those materials
separately; the controller neither downloads nor solves for replacements.

Supply the exact `pip-26.2.1-py3-none-any.whl` bootstrap artifact (1,816,632
bytes, SHA-256
`71138adf1f4ca900cdb7d289c21b7494329f2332b6d85f0e1c42108c0384ed3e`)
and a separately trusted, compatible Conda executable with its externally
recorded entrypoint digest. No Conda installation or its transitive bootstrap
dependencies are provided or authenticated by this component.

```bash
python3 -I -S -B /relocated/source-component/workflow/benchmark_tools/install_publication_base.py install \
  --receipt /relocated/source-component/workflow/benchmark_tools/results/reconstructed_base_fixture_20260927.json \
  --cache /absolute/acquired-base-archives \
  --conda /absolute/trusted-conda/bin/conda --conda-sha256 RETAINED_CONDA_ENTRYPOINT_SHA256 \
  --pip-wheel /absolute/acquired-wheels/pip-26.2.1-py3-none-any.whl \
  --output /absolute/fresh/base-preparation --acknowledge-historical-runtime
```

Use an output directory outside the immutable source component. Complete
input preflight precedes creation; execution has private HOME/cache/config,
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
