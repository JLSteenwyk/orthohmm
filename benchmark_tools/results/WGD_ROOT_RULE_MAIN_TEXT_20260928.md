# Root-Rule Diagnostic In Main Text

The main manuscript and current claim checklist now incorporate the completed
four-rule diagnostic from commit `55676f8a`. The text retains all three negative
alternative outcomes, unchanged native pair evidence and the distinction between
conditional software-rule effects and biological tree correctness. No parameter,
method default, reference cohort or benchmark score changed.

The [v8 PDF](publication_main_review_20260928_v8.pdf) and
[HTML review](publication_main_review_20260928_v8.html) supersede v7 for the
current main text. Historical exports remain unchanged. The main text now names
the controlled Threadripper panel explicitly rather than implying a separate
dedicated machine is required.

## Verification

- Fresh execution reproduces every retained diagnostic JSON field after
  serialization, including input checks and baseline partitions.
- Checks confirm four arms, six cases and seven families per arm; two identical
  alternative partitions, four eligible coverage losses under mapped-event and
  no recovery of any of the five focal homologs.
- All 66 focused tests pass across the renderer, diagnostic, trace,
  reconstruction, independent WGD score audit and scorer.
- The [render receipt](publication_main_render_20260928_v8.json) records 16
  citations, 17 local links and 16 tracked targets, with no Pandoc warnings.
- All five rasterized PDF pages were inspected. No clipping or incoherent overlap
  was observed; no word extends beyond a page boundary. The
  [layout/evidence receipt](publication_main_layout_20260928_v8.json) binds the
  revised manuscript, checklist, diagnostic and PDF.

This is evidence-backed manuscript integration, not new independent validation,
controlled timing or publication readiness. Full dependency/data-rights review,
unresolved QfO uncertainty and the remaining goal requirements stay open.
