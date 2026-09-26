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
