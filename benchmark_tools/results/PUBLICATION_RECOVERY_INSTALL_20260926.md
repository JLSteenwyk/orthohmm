# Fresh Recovery Environment Installation

A new private venv was installed entirely from a dedicated 11-wheel local
directory with `--no-index --require-hashes --only-binary=:all: --no-cache-dir`.
The [complete Python lock](publication_recovery_requirements_20260926.txt)
retains the prior audited wheels except for Leiden: version 0.12 is replaced
by the provider-digest-checked 0.11 wheel matching the recovery snapshot.
No existing environment, native executable or benchmark result was changed.

The [version-2 receipt](publication_recovery_install_20260926_v2.json) links
the preparation, pip installation report/log, full local byte audit, runtime,
smoke commands and output hashes. All 11 installed distribution identities
match the pip report, all report URLs point to the dedicated local wheelhouse,
and the wheelhouse inventory matches the report exactly. `pip check` passes.

The independent audit compares 3,052 installed package files byte-for-byte
against wheel contents and all 33 shipped scientific source files against
frozen Git revision `7f3a9e4`. All 15 Leiden distribution payload files match
the validated private snapshot. Isolated smoke probing requires imports and
dependency locations inside the new venv and successfully loads the three
OrthoHMM CPU libraries. Standard and high-sensitivity each return four groups
covering the same 38 genes exactly once; both partition SHA256 values are
`1115fd8193636510bbc8cc8462d1b874e2a50db662fc0a59d3552d811ffa0885`.

## Exceptions and Preserved Failure

The first byte audit rejected a NumPy `.pyc` rewritten during pip installation.
Its original auditor and failure message are retained. The corrected reader
excludes a bundled `__pycache__` member only when its corresponding wheel
source exists and passes the ordinary byte check. It does not prove bytecode
equivalence or weaken source/native-library matching.

Eleven top-level installed RECORD files, one regenerated NumPy bytecode file,
the relocated igraph C header and the numba command-line script are excluded
from byte-identity claims. The version-2 reader additionally checks 12 vendored
RECORD files previously excluded too broadly. The earlier receipt is retained
as superseded scope evidence, not the current audit. Generated entrypoints and
extra bytecode files are not an exhaustive installed-filesystem attestation.
Twenty-three focused wheel/install/Git-source tests pass, including rejection
of changed vendored metadata and missing or changed bytecode source.

## Executed Workflow

Artifacts are under `benchmarks/work/publication_recovery_install_20260926`.
An empty venv was created with `venv.EnvBuilder(with_pip=False)`. The existing
patched builder's pip installed into it using `--python`, with the lock above,
`--find-links .../wheels`, `--report .../install_report.json` and
`--log .../install.log`. The unchanged smoke CLI was then run:

```bash
python -m benchmark_tools.verify_cpu_wheel_install --root . \
  --python benchmarks/work/publication_recovery_install_20260926/venv/bin/python \
  --wheel benchmarks/work/publication_recovery_install_20260926/wheels/orthohmm-0.5.0-cp310-cp310-linux_x86_64.whl \
  --output /absolute/fresh/recovery-smoke
python -m benchmark_tools.audit_recovery_install --repo . \
  --directory benchmarks/work/publication_recovery_install_20260926 \
  --output /absolute/fresh/recovery-audit.json
```

The second command audits the retained successful smoke under the installation
directory. A separate smoke output is not silently substituted into that
evidence chain. Do not overwrite retained installations, logs or receipts.

This is a same-host Python-environment and small-fixture validation, not a
full OrthoBench reproduction in the newly installed environment. External
MAFFT/FastTree installation, canonical ordering integration, inferred-phylogeny
validation, transitive system dependencies, redistribution and cross-platform
testing remain separate requirements. The historical lock is not a claim that
all dependencies are current or secure. No publication readiness is asserted.
