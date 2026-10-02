# SwissTrees Descriptive Table Arithmetic

The four retained 2026-09-26 tables now have a standalone, standard-library
checker. It reproduces **208 rows across eight methods and 18 families** from
the admitted family sufficient counts and frozen bin assignments. The
[retained result](swiss_descriptive_table_arithmetic_20261002.json) checks
**984 score/difference cells**, including **744 numeric values and 240 missing
cells**. Each logical cell is checked in both JSON and TSV; these are not
1,968 independent endpoints. No retained score or bin assignment changes.

| Table | Rows | Score/Difference Cells | Missing Cells |
| --- | ---: | ---: | ---: |
| Entropy/length/explicit-description | 88 | 264 | 96 |
| Alignment identity | 32 | 192 | 48 |
| Historical fragment annotations | 56 | 336 | 48 |
| Mapped-tree duplication annotations | 32 | 192 | 48 |
| Total | 208 | 984 | 240 |

## Reproduce

From a checkout containing the seventeen committed payload files, run:

```bash
python -I -B benchmark_tools/check_swiss_descriptive_tables.py \
  --results benchmark_tools/results \
  --output /tmp/swiss-table-check.json
```

The output path must be new. Only Python's standard library is required;
the checker imports no project exporter, scorer or scientific library. It
uses fixed hashes for all four table manifests and verifies the recorded
counts, features, JSON/TSV rows and Markdown bytes at explicit current paths.
Embedded historical absolute paths remain provenance, not file-access targets.
All consumed payloads are rechecked after the final panel.

For each family, doubled directed counts with the retained pseudocount give
`P = (TP + 2)/(TP + FP + 4)` and `R = (TP + 2)/(TP + FN + 4)`.
The checker uses exact fractions, averages family precision/recall within
each bin, then calculates their harmonic mean. It does **not** average family
F1. Family statistics, full-family point estimates, all stratum estimates,
differences from full OrthoFinder, prediction semantics, family membership,
status, row/field inventories and missingness are checked. The fixed absolute
tolerance is 1e-12 on the proportion scale. Missing is not zero.

## Evidence And Limits

**49 local tests pass in 1.18s**, with zero failures/errors/skips. Tests include
malformed counts/admission/methods/bins, an explicit macro-versus-family-F1
counterexample, exact differences, nonfinite scores, fixed tolerance, metadata
types/semantics, missingness, changed bytes and a change after the last table
comparison. A standalone copied invocation succeeds and refuses overwrite.
Another fresh copied child rejects its original-checkout read canary, then
verifies all tables with zero original-path events; further subprocesses are
forbidden. All checked payload records and checker location are in that copy.
This is Python-level guarding, not OS containment or cross-host validation.
Earlier 44- and 47-case reports are retained, not added to the final count.

The seventeen payloads total 2,645,910 bytes and match committed source c78.
Source/report/test/log pins are in the
[validation receipt](swiss_descriptive_table_validation_20261002.json).
The four original exporters and their source/admission requirements remain
unchanged. Duplication export requires eleven retained source records;
fragment admission references 1,774 records totaling 607,096,715 bytes.
This derived-table checker does not pretend to revalidate those raw files or
fix their existing workstation-path CI failures. No annotation, alignment,
tree, biological inference, bootstrap, plotting or raw-QfO scoring is rerun.
Markdown is hash-checked, not independently regenerated or visually reviewed.

Separately, actual source-c78 macOS Python 3.11 fast CI confirms the previous
timeout panel: measurement 37, WGD eight and measurement replay ten all pass.
The full fast result is 13,765 pass, 29 fail, zero errors, 110 skip,
30 warnings, 362.00s. At 05:40:34 UTC all five source-c78 test jobs fail,
wheel/docs succeed. Only this fast log is inspected, and it predates the new
arithmetic checker. No sibling/new-checker/full-suite success is implied.

These are development-exposed descriptive associations, not causal effects,
independent biological validation or additional uncertainty estimates. Original
TreeFam sources, uncertainty for other QfO challenges, controlled resources,
raw-data/source rights, complete executable release and public deposition remain
open. Timing stays deferred; no contention poll, quiet-window question, DGX
access or action on unrelated jobs/services occurs. The publication goal remains
active. The existing review archive is not regenerated or expanded.
