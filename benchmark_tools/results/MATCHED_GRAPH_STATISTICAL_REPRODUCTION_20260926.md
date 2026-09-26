# Standalone Matched-Graph Statistical Reproduction

The standalone `benchmark_tools/reproduce_matched_graph_statistics.py` needs
only the retained `matched_graph_scores_20260926/results.json`, Python and
NumPy. It does not import OrthoHMM or other repository modules. It reconstructs
all per-dataset precision/recall/F1 values from TP/FP/FN, validates the complete
paired panel, and independently aggregates seed blocks and recomputes the
prespecified bootstrap effects/intervals. It checks all 24 metric contrasts,
eight adjusted F1 intervals and seed win/tie/loss counts within absolute
tolerance 1e-12.

```bash
python -I /relocated/reproduce_matched_graph_statistics.py \
  --results /relocated/results.json \
  --output /relocated/fresh_reproduction.json
```

A copy of the script and result JSON was run from `/tmp`, outside the
repository, with isolated Python and the separately installed NumPy 2.2.6.
The [reproduction receipt](matched_graph_statistics_relocated_20260926.json)
records the result, relocated files and file-access trace. A second traced
statistical invocation also passed. No native inference was rerun.
Exact quoted historical evidence paths from the result JSON were searched
in `strace -f -e trace=file`; none was accessed. This is a bounded file-access
check, not a claim that the Python/NumPy installation or OS runtime was moved.
The 17 focused tests include corrupted counts, effects, intervals, panel
membership and win counts, and all pass.

This closes a count-level reproduction dependency on the original repository
and simulation paths. It does not independently verify the truth, predictions,
graph inference, statistical sampling assumptions or biological generalization.
Only five reporting seed blocks were used; development exposure and limited
tail resolution remain. Python/NumPy are prerequisites, not archived here.
This is not a complete native workflow bundle, versioned publication release
or external deposit. The standalone script and exact result JSON are committed
alongside the report so the executable statistical check can be relocated.
