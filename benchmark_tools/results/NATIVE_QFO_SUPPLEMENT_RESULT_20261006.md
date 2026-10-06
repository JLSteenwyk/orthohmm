# Verified Native QfO Manuscript Supplement

The generated [seven-page PDF](native_qfo_manuscript_supplement_20261006_v3/supplement.pdf)
and [Markdown](native_qfo_manuscript_supplement_20261006_v3/supplement.md)
bring the two actual native P0/C0 cells, conditional SwissTrees intervals,
VGNC counts/changes, complete SwissTrees search/graph/reconciliation tracing,
unchanged sequence strata and functional-pair composition into one working
supplement. Historical selected-default comparator results remain separate.
Initial HMM search is on in both cells; this is not a total-HMM ablation.
No inference, scoring, scientific admission, bootstrap or failed-timing repair
is performed. Historical main text, PDFs and archive payloads are unchanged.

This is presentation progress, not full publication readiness. Five native
QfO identities are still absent from the supplied score snapshot. Conditional
SwissTrees F1 uncertainty includes zero, other-endpoint paired uncertainty
remains unresolved and R1 inference timing stays null/ineligible. Inferred
duplication labels, graph connectivity and descriptive feature associations
are not promoted into true evolutionary mechanisms or general superiority.
The supplement explicitly retains the unknown potentially tool-dependent
effects of shared-host contention; no quiet-window/dedicated-host gate returns.

## Generation And Checks

[Selection](native_qfo_supplement_selection_20261006.json) pins 20 existing
reports/manifests/readbacks. The generator checks their current sources and
cross-report bindings, the complete seven-identity/two-admission cohort,
null missing scores, failed R1 timing, unchanged 100,000 draws/seed20260922/
42-endpoint adjustment, diagnostic cohort agreement and observed-artifact
admission limits. These are direct checks, not repeated transitive scientific
admissions or raw-data scans. The selection is unchanged through all attempts.

Eight generated TSVs retain 51 rows, including zero transitions and empty bins.
The manuscript displays 34 rows across eight tables; its SwissTrees transition
display omits zero cells and its sequence display uses the original entropy/
relative-length bins. Complete rows remain in the TSVs and source reports.
[Assembly receipt](native_qfo_manuscript_supplement_20261006_v3/assembly.json)
records input, source and output identities, without a goal-complete flag.

The first actual invocation at `9888cd1e` refuses a fixture/production schema
type mismatch: `new_bootstrap_draws` is an integer count zero, not boolean false.
No output directory or scientific result is produced. Preserve the
[initial refusal](NATIVE_QFO_SUPPLEMENT_INITIAL_REFUSAL_20261006.md), original
source revision and failed GNU-time receipt. Fix `a780d653` requires integer
zero and adds boolean/float/negative-count refusals, without changing intervals.
Initial 30-test and corrected 33-test receipts remain local.

Version2 then generates and prints six pages successfully, but agent visual
inspection finds repeated implicit Pandoc figure captions. Preserve version2
Markdown/HTML/PDF, tables, receipts and previews as presentation history, not
a failed scientific run. Fix `3f457e78` uses a short image alt and inline caption
and adds a no-duplicate-caption test. All eight version2/version3 TSVs are
byte-identical. No scientific computation is replayed for either layout.

The [independent presentation check](native_qfo_manuscript_supplement_20261006_v3/presentation_readback.json)
imports no generator. Structured Pandoc AST parsing reproduces all eight
displayed table bodies from CSV-read TSVs, including declared formatting and
display filters; all 34 rows also occur in printed PDF text. Both embedded
figures decode to exactly the original PNG dimensions/channels/pixel hashes,
not just nonblank rectangles. This verifies presentation correspondence, not
an independent raw scientific re-admission. Both original scientific figure
readbacks remain bound, and no figure is redrawn.

[HTML asset receipt](native_qfo_manuscript_supplement_20261006_v3/html_assets.json)
checks 33 local occurrences/24 unique already tracked targets. Chrome prints
with its default sandbox, not `--no-sandbox`; its launcher identity is recorded,
not a complete browser-runtime closure. The copied PDF matches the original
[print receipt](native_qfo_manuscript_supplement_20261006_v3/print_receipt.json).
[Layout review](native_qfo_manuscript_supplement_20261006_v3/layout_review/report.json)
has no bounds violations and renders all seven pages. Earlier viewed pages1-3
are byte-identical to final renders; final pages4-7 were separately viewed.
Tables/figures/captions are readable and no incoherent overlap was observed.
This is agent visual inspection, not a human-review certification; automatic
receipts correctly keep visual-review flags false.

Final joined tests:143 passed in4.48s, no failures/errors/skips. This includes
33 generator and six presentation-check cases plus existing render/print/PDF
and VGNC contracts. Retain earlier137/143-test receipts and final local XML
`benchmarks/results/native_qfo_supplement_joined_tests_20261006_v3_final.xml`.
The final generator runs in original Python3.10.13; rendering/checks use the
separate Python3.12/PyMuPDF1.27.2.3 environment. No scientific environment is
installed into or modified. GNU-time final assembly0.24s/64,576KiB/zero swaps
is shared-host postprocessing, not native inference timing or tool-speed evidence.
About765GiB RAM was available, with near-full host swap retained as a caveat.

## Identities And Reproduction

Selection SHA256 `241f9753f9a360fcc0ae84ff14c152f518982104f020d25a0c26f9eaa3237c4a`.
Generator SHA256 `999a5a5ceb0a19b4467c5c20cb059ee4eaddd4c954a89bfa1ff2a8fb26750e07`.
Presentation checker SHA256 `3fae7001026b28b4605e5b0d0c9f2465e2932597cff3c6a54f638236dbfe7d61`.
Markdown11,851bytes SHA256 `10e09bba834aa94679661342e112660aa355c3a94d3660b0498a3dec2ab44200`.
PDF563,905bytes SHA256 `c5e82c28e4dd04e47c4666560171b05494dfac24969df011e6a8f0965ddafc04`.
Assembly22,153bytes SHA256 `1dcad56e356bd7f25304b10e5545adbc5117b4badda250c0272000688f37167d`.
Presentation readback3,483bytes SHA256 `2bc8e5b90b7a994466cde0f526ac11d4f63f4cdf30c4fda02fd879a18ee7b8e8`.

Generate to a fresh repository-local output directory:

```bash
env -u PYTHONHOME -u PYTHONPATH -u PYTHONUSERBASE -u LD_PRELOAD \
  -u LD_LIBRARY_PATH -u LD_AUDIT PYTHONNOUSERSITE=1 \
  PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  benchmark_tools/assemble_native_qfo_supplement.py --root . \
  --selection benchmark_tools/results/native_qfo_supplement_selection_20261006.json \
  --selection-sha256 241f9753f9a360fcc0ae84ff14c152f518982104f020d25a0c26f9eaa3237c4a \
  --output benchmark_tools/results/native_qfo_supplement_replay
```

Use existing `render_manuscript_review`, `print_manuscript_review` and
`review_manuscript_pdf` CLIs with fresh paths; the retained receipts contain
actual commands. `check_native_qfo_supplement_presentation --help` gives the
assembly/assets/print/PDF/output arguments for checking a new rendering.
Local absolute source/evidence paths are still required. This is not portable
full-study restoration, a refreshed final archive or external deposition.

Original22444 remained RUNNING8:29:35;22445/22450/22451/22452 dependency-pending.
No unfinished output was read, original handle restarted, native successor
released or unrelated workload changed. Remaining native cells, matched search,
wider uncertainty/generalization/error strata, original TreeFam/provenance and
full manuscript/archive/release work remain required. Full goal active.
