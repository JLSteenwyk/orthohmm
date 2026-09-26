# Relocated Forced-Candidate Figure Reproduction

Exported committed sources and the admitted count summary at revision
`002340e46943fc38593426b494e952061bbb36c2` to an empty directory outside
the repository. Ran the original plotting module there using Python's
isolated mode, adding only the relocated export to its import path.
All three formats were generated. The PNG is byte-identical to the retained
figure (`cmp` exit zero), SHA256
`0612cabb4522f653fcc9be672447f73e54666eadf8a49de984fe93af49636a4b`.

The local archive is `benchmarks/work/forced_candidate_figure_sources_20260926.tar`,
SHA256 `5226092bd7e1760917d660c67b677fffb83dda5f449bb7ff994f2942c59d8f22`.
It contains nine Python modules, the audit JSON and repository license,
not proteomes, raw searches or reference trees. Recreate it from Git:

```bash
git archive --format=tar --output=/tmp/forced-figure-sources.tar \
  002340e46943fc38593426b494e952061bbb36c2 \
  LICENSE.md benchmark_tools/__init__.py \
  benchmark_tools/benchmark_production.py \
  benchmark_tools/build_publication_runtime.py \
  benchmark_tools/plot_forced_candidate_diagnostic.py \
  benchmark_tools/prepare_ob_candidate_neighborhood.py \
  benchmark_tools/run_simulation_generation.py \
  benchmark_tools/run_simulation_methods.py \
  benchmark_tools/validate_profile_runtime.py \
  benchmark_tools/verify_simulation_histories.py \
  benchmark_tools/results/ob_forced_candidates_audit_20260926.json
```

Extract into a fresh directory and run there with Python 3.10.13 and
Matplotlib 3.10.8 (the tested existing environment):

```bash
python -I -B -c 'import pathlib,runpy,sys; sys.path.insert(0,str(pathlib.Path.cwd())); sys.argv=["plot", "--audit", "benchmark_tools/results/ob_forced_candidates_audit_20260926.json", "--output", "rerendered"]; runpy.run_module("benchmark_tools.plot_forced_candidate_diagnostic",run_name="__main__")'
```

The audit's historical paths remain provenance and are not opened by plotting.
This verifies relocated plotting and count-marginal checks, not a rerun of
HMM search, independent reconstruction of the raw counts, or installation
from an empty environment. PDF/SVG were generated but are not claimed
byte-identical because format metadata can differ. The complete publication
workflow, dependency lock, rights review and external archival deposit remain
separate requirements.

## Clean Environment Follow-Up

A fresh venv, using Python 3.10.13 with `include-system-site-packages = false`,
also regenerated the identical PNG and passed `pip check`. Eleven dependencies
were installed from PyPI wheels with explicit versions. The resulting
[Linux/Python 3.10 requirements](forced_candidate_plot_linux_py310_requirements_20260926.txt)
pin the SHA256 of each tested wheel. Install them into a new venv with:

```bash
python -m pip install --index-url https://pypi.org/simple --only-binary=:all: \
  --require-hashes -r forced_candidate_plot_linux_py310_requirements_20260926.txt
```

Use that venv's Python for the isolated plotting command above. This follow-up
adds a clean plotting-install test; it still does not reproduce HMM inference
or independently reconstruct raw benchmark counts. The wheel lock is specific
to the tested CPython/Linux x86_64 platform, not a cross-platform lock.

## Pillow Security Update

The [read-only GitHub alert snapshot](dependency_alerts_release_review_20260926.json)
identified 13 Pillow advisories affecting the initially reproduced 12.2.0
wheel. All list 12.3.0 as the first patched version. Updated the plotting
lock to Pillow 12.3.0 and its tested wheel hash, installed that wheel in the
isolated venv, and regenerated all three formats. PNG comparison against
the original figure again passed byte-for-byte. This patch does not change
the scientific benchmark runtime or stored evidence. The initial historical
installation remains documented, but is no longer the recommended lock.
The snapshot also lists 11 pip/setuptools alerts in the separate historical
CPU wheel manifest; this Pillow update does not resolve those alerts or
constitute a full installed-environment security audit.
