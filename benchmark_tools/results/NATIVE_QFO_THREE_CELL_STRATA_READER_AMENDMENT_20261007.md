# Three-Cell Readback Path-Type Amendment

The exporter committed/pushed at a5952bbc completed once with exit0 and
54family/60score/40contrast rows. Its report55ed72ce and all outputs are
preserved. The first independent reader failed with exit1 at TSV path handling:
JSON output paths are strings but the helper called path.open(). No readback
output was written. Retain the original sourcef59bd8cc and the observed
failure receipt; do not edit that failed reader or rewrite the export.

Prospective correction: new readback_native_qfo_three_cell_strata_v2.py uses
Path(path).open() in the TSV helper. All input checks, exact rational arithmetic,
inventories, table checks and scope requirements remain unchanged. Add fixture
tests with actual JSON string paths and a full temporary export/readback
workflow, not just individual helpers. After these tests, commit/push the
new reader, tests, this amendment and the original failure receipt before
one fresh v2 readback of the unchanged report. This is a diagnosed code
repair and versioned validation attempt, not an automatic failed benchmark
retry or a rerun of inference/scoring/subgroup calculations.

The v2 reader binds this amendment and failed source/receipt in its checked
inputs. Preserve any further actual failure without overwrite. No new bins,
defaults, evidence admission, scientific timings, intervals or biological
claims. Existing23902/23910 remain unchanged and goal remains incomplete.
