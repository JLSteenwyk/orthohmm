# Relocated Descriptive Figure Reproduction

Exported six committed files from revision
`0e072fce5defeff363c3d0a03176e92f514f2c2c` into
`/tmp/orthohmm-descriptive-features-20260923`, outside the checkout: the
renderer, its standard-library provenance helper, two score tables, expected
endpoint table and license. All six exported byte strings were independently
compared with `git show` at that commit after execution.

Executed the renderer with Python isolated mode (`-I -B`) in the existing
`swiss_analysis_env_20260917` virtual environment: Python 3.10.13, NumPy 2.2.6,
Matplotlib 3.10.8. The clean run explicitly unsets LD_PRELOAD, LD_LIBRARY_PATH,
LD_AUDIT, PYTHONPATH and PYTHONHOME. Only the exported directory is explicitly
added to the import path; no full repository package is supplied. The virtual
environment prefix differs from its base interpreter prefix.

The [verification receipt](swiss_descriptive_feature_reproduction_20260923.json)
records exact agreement for all 144 endpoint rows, including 18 unavailable
values. Endpoint TSV SHA256 is
`1e8cabb0bbec77396b04c36fc536cc7562c29f1a36cf0950e79b6acb188c9f38`.
Both 2700-by-1620 PNGs are nonblank and their RGB pixel arrays exactly match
the committed originals. PDF/SVG file-byte equivalence is not claimed;
renderer metadata may differ. Their generated hashes are retained.

An earlier exploratory run in the same isolated Python environment inherited
a CUDA LD_LIBRARY_PATH. Its output remains under `reproduced`; the clean
follow-up is under `reproduced_clean`. Only the clean follow-up is retained
as the receipt's generated-artifact panel. No earlier result was overwritten.

## Reproduction Commands

Use a fresh export directory and set `PYTHON` to an isolated environment with
the versions above. The six-file export can be recreated with:

```bash
git archive 0e072fce5defeff363c3d0a03176e92f514f2c2c \
  LICENSE.md \
  benchmark_tools/plot_swiss_descriptive_features.py \
  benchmark_tools/prepare_ob_candidate_neighborhood.py \
  benchmark_tools/results/swiss_identity_strata_20260923/scores.tsv \
  benchmark_tools/results/swiss_fragment_strata_20260923/scores.tsv \
  benchmark_tools/results/swiss_descriptive_feature_figures_20260923/endpoints.tsv \
  | tar -x -C "$EXPORT"
env -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT -u PYTHONPATH -u PYTHONHOME \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0 \
  "$PYTHON" -I -B -c '
import sys, runpy
from pathlib import Path
root = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(root))
sys.argv = ["plot", "--root", str(root), "--output", str(root / "reproduced_clean")]
runpy.run_module("benchmark_tools.plot_swiss_descriptive_features", run_name="__main__")
' "$EXPORT"
cmp "$EXPORT/benchmark_tools/results/swiss_descriptive_feature_figures_20260923/endpoints.tsv" \
    "$EXPORT/reproduced_clean/endpoints.tsv"
```

This is a same-host visualization reproduction from frozen descriptive
tables, not cross-platform validation, fresh annotation acquisition,
statistical estimation, native inference, raw scoring, a complete archival
bundle or publication readiness. Raw input rights are not resolved by this
bounded export.
