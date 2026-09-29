# Private Timing Runtime Candidate

Update: the [packaging alignment](THREADRIPPER_PACKAGING_ALIGNMENT_20260928.md)
subsequently replaced only private packaging 26.0 with hash-proven historical
26.1 payloads. The comparison is now 709/710 identical files. The original
installation/lock below is historical; use the linked v2 lock for that revision.
The PyYAML build difference and native timing integration remain open.

Created a separate candidate at
`benchmarks/work/threadripper_private_runtime_20260928/venv`, using the already
reconstructed CPython 3.10.13 at
`/tmp/orthohmm-base-reconstruction-20260927/python-runtime`.
No shared package was installed, removed or changed. The venv uses copied
executables and no system site-packages, but still depends on that private
base prefix and shared OS libraries; it is not a self-contained container.

The [installation receipt](threadripper_private_runtime_20260928.json)
records 30 binary wheels (164,266,582 bytes), versions, digests, commands and
logs. Twenty-nine packages use version labels from the frozen timing baseline,
including Numba 0.65.0, llvmlite 0.47.0, DendroPy 5.0.8, igraph 1.0.0,
leidenalg 0.11.0 and pyparsing 3.2.1. Pip 26.2.1 is explicitly an installer
selection. The separate recovery runtime's newer Numba/DendroPy versions
were not substituted. Historical pins are not security recommendations.

Wheels were downloaded from the explicit public PyPI index without source
builds. The generated [hash lock](threadripper_private_requirements_20260928.txt)
was then installed offline with required hashes and no dependency resolution.
All five stages completed: download, venv creation, installation, `pip check`
and declared native imports. Installed distribution inventory matches the
selected labels exactly; `pip check` reports no broken requirements.

## Import Comparison

The native import probe explicitly loads the frozen source directory under
isolated Python mode. No imported module or Python search path contains the
shared `/home/bizon/anaconda3` prefix. This is an observed Python lookup check,
not a file-access sandbox or native-library closure audit.

The [byte comparison](threadripper_private_import_comparison_20260928.json)
rechecks all 40 frozen core-file records and compares 710 imported site-package
or scientific-core module files with the historical lookup: **707 match**.
All ten pyparsing files now match their historical hashes. Six unrelated
editable-startup finder modules are intentionally absent. Two additional
module names are explicit in the result, including `packaging._structures`.

Three file differences prevent treating this as an identical runtime:

- `packaging` and `packaging.version`: the baseline distribution inventory
  lists 26.0, but the historical imported `packaging/__init__.py` declares
  26.1. Its current shared copy still matches the historical hash, so this
  discrepancy predates the present candidate. The candidate follows the
  metadata label 26.0 and therefore does not reproduce those module bytes.
- `yaml._yaml`: the PyYAML 6.0.1 wheel extension differs from the historical
  extension; matching version labels do not establish binary identity.

The candidate is retained unchanged for review, not promoted to the timing
baseline. Resolve the packaging source/metadata mismatch, review the YAML
build and native linkage, and validate lookup plus native output parity under
an explicit deployment amendment. Do not infer prediction or runtime effects
from these differences alone. No fixture or production timing run was launched.

25 focused tests pass for staging, comparison and preceding drift tools.
Cases cover version selection, wrong pins, duplicate/missing/wrong-version
wheels, ambiguous metadata, changed/omitted modules and no-overwrite behavior.
The wheelhouse remains local; no binary redistribution clearance, complete
release or timing eligibility is asserted. Final accounting, environmental
validation and a quiet window remain separate requirements.
