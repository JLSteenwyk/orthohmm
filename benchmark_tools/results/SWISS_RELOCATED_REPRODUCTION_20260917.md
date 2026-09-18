# Relocated SwissTrees Reproduction

Current reproduction uses the patched Pillow12.3.0 analysis lock; see
[security remediation and exact rerun](SWISS_ANALYSIS_SECURITY_20260917.md).
The initial environment and hashes below are retained historical evidence,
not instructions to reinstall the vulnerable Pillow12.2.0 pin. Commands use
the current lock at checkout.

The completed SwissTrees comparator statistics and figure workflow were rerun
from a clean, explicit source export of commit `fd5f1ac`, outside the project
worktree, in a fresh Python 3.10.13 virtual environment. This is a statistical
workflow reproduction, not native inference or raw-QfO benchmark reproduction.

The export contains nine committed analysis modules, four committed result/
protocol files and the repository license. It does not copy the dirty worktree,
raw QfO installation, inference outputs or unrelated experiments. The runner
uses isolated Python mode, disables user-site imports, removes Python/loader
overrides, uses a fresh plotting configuration directory and fixes numerical
library threads to one. Embedded historical paths remain provenance only.

## Environment

`benchmark_tools/swiss_analysis_requirements.txt` pins the 11 installed packages;
`swiss_analysis_requirements.lock` includes distribution hashes. It was generated
with uv 0.12.15 and installed using hash verification into a new virtual
environment. `uv pip check` reports all 11 packages compatible. This environment
is separate from every frozen inference environment. It is not a security audit
or cross-platform compatibility claim.

Lock SHA256: 12e6f4dc08c0566b269c7f0c30f70f0b1e2bc55e220424e1262c5702146cbfff.

## Reproduce

From the repository root, with uv 0.12.15 and Python 3.10 available, choose new
environment, export and report paths:

```bash
uv venv --python python3.10 /tmp/orthohmm-swiss-env
uv pip sync --python /tmp/orthohmm-swiss-env/bin/python --require-hashes benchmark_tools/swiss_analysis_requirements.lock
uv pip check --python /tmp/orthohmm-swiss-env/bin/python
python benchmark_tools/reproduce_swiss_comparators.py --revision fd5f1ac --python /tmp/orthohmm-swiss-env/bin/python --output /tmp/orthohmm-swiss-reproduction --report /tmp/orthohmm-swiss-reproduction.json
```

The workflow refuses existing output directories and reports. The source
revision must exist locally. Installation can require package-index access;
dependencies are not bundled for offline installation. The Python interpreter
itself and system libraries are recorded/tested here, not shipped in a container.

## Observed Results

- All scientific JSON content matches exactly, including all 100,000-draw
  intervals, point estimates, family contrasts, descriptive counts, controls,
  NumPy version and limitations. Only relocated provenance paths may differ;
  corresponding input/helper/source bytes and hashes must remain identical.
- The generated Markdown table is byte-identical to the committed table.
- The pinned-report plotter runs successfully after statistical equivalence is
  established, producing PDF, SVG and PNG. The relocated PNG was visually
  inspected and all panels, labels and intervals render without overlap.
- No bitwise PDF/SVG identity is asserted; renderer metadata may vary.
- All 21 focused reproduction, figure and bootstrap tests pass. Mutation tests
  reject changed intervals, missing family records, version changes, input or
  helper identities and unexpected result fields.

Machine-readable evidence: `swiss_relocated_reproduction_20260917.json`, SHA256
2df2a6226fd4a4c56d0ce8076e87d8677d6c34a3c41a9f94d1b4ed9b38ad31ca.
It records the resolved commit, exported files, environment, executed commands,
logs and generated artifacts. The initial local package installation could not
hardlink across filesystems and fell back to copying; installation and dependency
checks succeeded. No analysis result was substituted or altered.

This closes a bounded relocation check for one completed analysis. The broader
publication package still needs other executable workflows, raw-data acquisition
and scoring reproduction, external-data licensing review, complete inference
environments, the final versioned release and archival deposition. Exporting the
repository license does not resolve third-party data redistribution rights.
