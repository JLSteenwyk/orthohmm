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

Dataset/reference files, test samples, native predictions, scoring outputs,
result receipts/plans, figures, manuscript assets, third-party source archives,
dependency wheels, binaries and runtime/OS images are excluded. This avoids
newly redistributing those materials; it does not resolve every source-file
attribution or clear an arbitrary future release. Both project license copies
are retained. Do not interpret their presence as an override of third-party
notices or terms.

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
