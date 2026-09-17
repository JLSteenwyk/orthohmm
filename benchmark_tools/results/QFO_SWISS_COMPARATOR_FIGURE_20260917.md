# SwissTrees Comparator Figure

[PDF](figures_qfo_swiss_comparators_20260917/swiss_comparator_intervals.pdf),
[SVG](figures_qfo_swiss_comparators_20260917/swiss_comparator_intervals.svg),
[PNG](figures_qfo_swiss_comparators_20260917/swiss_comparator_intervals.png).

The three panels display F1, precision and recall differences for all eight
prespecified contrasts. Points show observed differences; thick lines show
nominal 95% paired percentile intervals; thin lines show Bonferroni-adjusted
intervals across 24 contrast/metric endpoints. All plotted values are converted
from raw 0-to-1 differences to percentage points by multiplication by 100.
Identical axes are used across panels and zero is explicitly marked.

The first seven rows compare each candidate with full OrthoFinder 3.1.5. The
separated last row compares phylogenetic and high-sensitivity OrthoHMM. The
MCL-checkpoint and supplied-tree FastOMA diagnostics are explicitly labeled.
The final row is not a pure reconciliation ablation.

Input: `qfo_swiss_comparator_intervals_20260917.json`, SHA256
6317d0274142354dbc0745dc3c75a2532d0d39dedeba8f30b9aca8c6fee21005.
The artifact manifest records source/input/output hashes and plotting versions.
No new statistics or inference runs are introduced by plotting.

Reproduce from the repository root with a fresh output directory:

```bash
python benchmark_tools/plot_qfo_swiss_comparators.py --results benchmark_tools/results/qfo_swiss_comparator_intervals_20260917.json --output /tmp/orthohmm-swiss-comparator-figure
python -m pytest tests/unit/test_plot_qfo_swiss_comparators.py -q
```

Validation compares every plotted point and interval endpoint with source values,
checks retained contrast order and resampling controls, rejects nonfinite,
unordered or clipped intervals, and checks text bounding boxes against the
canvas. The generated PNG was also visually inspected: all panels, contrasts,
labels and caveats are visible without overlap. Five focused figure tests pass.

These approximate intervals describe 18 curated, development-exposed families.
Shared history and merged predictions can correlate families; the figure is
not independent validation or evidence about other QfO metrics or the secondary
six-metric summary. See [full results](QFO_SWISS_COMPARATOR_INTERVALS_20260917.md)
and the [frozen protocol](QFO_SWISS_COMPARATOR_UNCERTAINTY_PROTOCOL_20260917.md).
