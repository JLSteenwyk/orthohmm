# All-Method OrthoBench Strata Results

Protocol and exporter were committed and pushed as `0a28a042` before joining
the five additional methods to the frozen bins. The export contains all 112
rows: fourteen strata for eight methods. All 42 original three-method rows
agree within 1e-10 percentage points. Sixteen rows represent two empty bins
and remain missing, not zero.

The [complete F1 table](ob_complete_strata_20260928/TABLE.md) includes every
stratum. [Machine-readable scores](ob_complete_strata_20260928/scores.tsv)
also retain precision and recall. The [report](ob_complete_strata_20260928/report.json)
preserves family membership, weighted sufficient statistics and provenance.
A separate standard-library checker recomputed the weighted counts and
ratios with exact rational arithmetic, checked the TSV and rehashed all 185
source-file records. Its [receipt](ob_complete_strata_20260928/crosscheck.json)
binds the export and checker bytes.

These are development-exposed descriptive results, not new inference or
independent validation. No confidence intervals, rankings, significance tests
or selected winning strata were added. Existing three-method intervals and
the 84-endpoint correction do not apply to added comparators. Small bins,
including the one-family concentrated-composition bin, must not support broad
claims. Relative length, copy number and composition remain descriptors,
not validated fragmentation, duplication-history or domain annotations.

From the repository root, using an unused export directory:

```sh
python -B -m benchmark_tools.export_ob_complete_strata --root . --output benchmark_tools/results/ob_complete_strata_reproduction --protocol-sha256 3534b2839c8907c7affd445c0bf23e572b6d39654c15f130691192f1c9e52184
python -B -m benchmark_tools.check_ob_complete_strata --base benchmark_tools/results --export benchmark_tools/results/ob_complete_strata_reproduction
```

The checker must run without Python's `-O` option, which disables assertions.
Retained source records contain local absolute paths; this is an in-place
reproduction, not a standalone portable raw-data bundle. No native inference,
new official-scorer invocation or historical input-consumption proof is
provided. No scientific defaults or timing measurements changed.
