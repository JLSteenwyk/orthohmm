# Complete SwissTrees Portable Arithmetic Reproduction

All 24 endpoints of the completed eight-method comparator analysis reproduced
within 1e-12 outside the repository, in Python isolated mode. The
[receipt](swiss_complete_portable_reproduction_20260926.json) binds the
retained count report and verifier hashes. This is statistical arithmetic
reproduction, not native inference, raw-count admission or a full release.

## Executed Export

Exported exactly three tracked files from commit
`6fa3473a1b038a07802cd1a3014bee8adc1f3da1`: `LICENSE.md`,
`benchmark_tools/reproduce_corrected_swiss_comparison.py`, and
`benchmark_tools/results/qfo_recovered_swiss_uncertainty_22178.json`.
Extracted the archive to `/tmp/orthohmm-complete-swiss-ufXcEi`.
The first archive attempt used the nonexistent path `LICENSE` and failed;
the corrected export uses `LICENSE.md` and completed successfully.

The retained local archive is
`benchmarks/work/complete_swiss_arithmetic_sources_20260926.tar`, SHA256
`558fc5bc6b7d193de905dad4887eb8568ae16adf01d2bc76a24f0560e448eb83`.
It can be regenerated with `git archive` from the revision and three paths
above; no raw sequences or native trees are required for this arithmetic check.

Executed from the extracted directory using the already isolated plotting
environment's Python 3.10.13 and NumPy 2.2.6:

```bash
/tmp/orthohmm-forced-figure-OQYDVY/clean-env/bin/python -I -B \
  benchmark_tools/reproduce_corrected_swiss_comparison.py \
  --results benchmark_tools/results/qfo_recovered_swiss_uncertainty_22178.json \
  --results-sha256 a599749d66433ec211ed1ab0e3a6a2a75eb99cb0dc57abc47c2d267930cbbe34 \
  --output reproduced.json --retained-counts-only
```

The [existing wheel lock](forced_candidate_plot_linux_py310_requirements_20260926.txt)
pins NumPy and the plotting environment; only NumPy is needed by this
standalone verifier. This run reused that clean environment, rather than
claiming a second newly installed environment.

## File-Access Check

Repeated the same command under `strace -f -e trace=%file`, with a fresh
`traced_reproduced.json` output. The two receipts are byte-identical.
The trace contains no occurrence of the original repository's absolute path,
`qfo_benchmark`, or `benchmarks/work`. The exported script only invokes its
retained-count branch; historical provenance paths are not dereferenced.
The local trace is
`benchmarks/work/swiss_complete_portable_file_access_20260926.log`, SHA256
`3209c2d62e52a1ad78c7e58f325f468617806fee2838474546af739adc811e36`.

Both executions completed with exit 0. This checks relocation on the same
host and the same NumPy RNG/quantile implementation. It does not establish
cross-platform equivalence, independent reference correctness, redistribution
clearance, native workflow portability or publication readiness.
