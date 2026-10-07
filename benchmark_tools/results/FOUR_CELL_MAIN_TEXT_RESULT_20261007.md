# Four Cell Main Manuscript Integration

The new [v4 main manuscript](PUBLICATION_MAIN_TEXT_20261007_v4.md) integrates
the four admitted native QfO configurations, three supported conditional
SwissTrees contrasts and independently checked profile-pair localization.
All-tool comparisons, frozen methods, prior P0 error strata and citations
retain their earlier scientific scope. Three missing score cells and 11
missing contrasts remain unavailable; failed timing is not repaired.
All three adjusted SwissTrees F1 intervals include zero. This is neither
independent confirmation nor evidence for promoting a default.

## Generation And Content Checks

Source, tests and prospective protocol were committed/pushed as `e7c1145d`
before one selected generation. Clean scientific Python 3.10 execution:

```bash
env -u PYTHONPATH -u PYTHONHOME -u PYTHONUSERBASE \
  -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /usr/bin/time -f 'elapsed_s=%e max_rss_kib=%M exit=%x' \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  -m benchmark_tools.prepare_four_cell_main_text \
  --evidence-directory benchmark_tools/results \
  --output benchmark_tools/results/PUBLICATION_MAIN_TEXT_20261007_v4.md \
  --receipt benchmark_tools/results/publication_four_cell_main_generation_20261007_v1.json
```

Actual exit 0, elapsed 0.05 seconds and maximum RSS 12,288 KiB. Stdout:

```json
{"native_cells": ["p0_c0_r0", "p0_c0_r1", "p0_c1_r0", "p1_c0_r1"], "endpoint_rows": 6, "coverage_rows": 4, "conditional_interval_rows": 9}
```

The [generation receipt](publication_four_cell_main_generation_20261007_v1.json)
has SHA256 `2bc797440045bbe03c3c4ec011403abfa409c4272cccdd4f501cc305bd50ea30`.
It binds selected source `7c09719e4de2bd8182fd5a218ee4249b175ead6ad39b6aceddefca8362a78789`
and generated manuscript `0b7012e3dd46bd3b028857bdb3cbe3a9e9048936b90f3df9e4833cb342d5705a`
(74,287 bytes). Selected generator and manuscript are now immutable.

A separate readback of the actual generated Markdown used Pandoc's table AST
and direct snapshot/binding access, without invoking generation. All 24 scores,
four coverage/count rows and nine interval/effect/win-tie-loss rows match.
Removing the replaced native section and reversing only the three documented
ancillary clarifications restores the entire original scientific body exactly.
The actual v3-to-v4 diff was inspected. Original v3 and rendering sources have
empty scoped diffs. The 29 focused tests passed again in 0.78 seconds before
generation; the earlier XML remains unchanged.

## Review Scope

At this generation checkpoint no new HTML, printed PDF, visual review or
archive inclusion is established. Subsequent review evidence must identify
this exact v4 source and use fresh namespaces, not relabel an older review.
The previous v3 review, historical figures, source snapshots and archives
remain unchanged. No new raw inference, scoring, bootstrap, admission or
scientific default selection was performed. The full publication goal remains
incomplete, including actual pending native outcomes and disclosed limitations.
