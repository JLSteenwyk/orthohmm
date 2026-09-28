# Reporting Collector Native Outcomes

Jobs 22367-22369 completed 0:0 sequentially without restart or retry. The
[outcome receipt](threadripper_reporting_fixtures_outcomes_20260928.json)
pins terminal controller observations, driver/result records, independent
resource/output audits, release guards, reporting records and output parity.
All pre/post runtime trees match the v4 binding. Native lookup signatures
reproduce the pinned baseline for both interpreters in all six checks.
Release-budget decisions reproduce from their retained controller text.

| Job | Method | Reporting seconds | Job peak through reporting, bytes | Maximum foreign cores |
| --- | --- | ---: | ---: | ---: |
| 22367 | OrthoHMM high sensitivity | 0.257785149 | 726,667,264 | 71.3577 |
| 22368 | OrthoHMM satellite_v2 | 0.288694906 | 726,663,168 | 74.5699 |
| 22369 | OrthoFinder full | 0.257567402 | 360,677,376 | 74.6587 |

These are 16-gene integration fixtures on a busy host, not comparative
performance results. The peak values include preparation/verification and
observer work through the reporting read, but exclude later validation and
teardown. Every post-report peak equals its post-native peak in this small
fixture; that does not establish negligible reporting memory on large runs.
No native/observer time or memory was subtracted, and no controlled timing
was admitted. Reporting records independently validate under schema v4.

## Exact Output Parity

Input-species mappings and group partitions match previous jobs 22363, 22364
and 22365, respectively. All contain 16 genes in three groups. Satellite_v2
root HOGs match too. Its 36 canonical native ortholog pairs and OrthoFinder's
36 canonical native pairs match their corresponding prior runs exactly.
Native format/duplicate checks are separate from canonical-set fingerprints
and also passed. This demonstrates behavior preservation on the fixture,
not genome-scale equivalence or biological accuracy.

Seven additional [corrupted-record controls](threadripper_reporting_rejection_controls_20260928.json)
against the completed high-sensitivity fixture were rejected: failed status,
wrong job, wrong scope, wrong reporting duration, an early memory read,
a decreased/invalid memory peak and a missing memory field. These changed
JSON reads in memory only; raw records remained untouched. This is bounded
negative testing, not exhaustive validation. The focused renderer, reporting,
runtime, driver and output-audit test selection also passed (49 tests).

## Manuscript Export

The refreshed [five-page main draft](publication_main_review_20260928_v4.pdf)
and [HTML](publication_main_review_20260928_v4.html) now include the already
written FAS sampling limitation: scored fractions span 0.0067%-58.46% and
samples reuse proteins across pairs. Scientific source text was unchanged.
The [render receipt](publication_main_render_20260928_v4.json) records 13 local
link occurrences, 12 unique tracked targets, 16 resolved citation identifiers
and no pandoc warnings. The [layout receipt](publication_main_layout_20260928_v4.json)
records source hashes, PDF/HTML identities and all five rasterized pages.

All five pages were visually inspected; text and bibliography are readable
without visible clipping or overlap. Automated word-bound checks found no
page-bound violations at a one-point tolerance. The FastOMA bibliography
entry spans pages 4-5; journal-specific pagination is still pending. This is
a working review export, not an archival bundle or scientific validation.
The previous v3 export remains unchanged as a historical snapshot.

## Next Work

All three diagnostic jobs are terminal. No production timing identity was
started; the 27-run panel still requires integrated orchestration, a verified
quiet window and full-scale observer validation. Complete-job resource
accounting beyond report generation is also unfinished. Independent-family
confirmation, remaining uncertainty and the broader release/publication
requirements remain open. No unrelated work or DGX service was changed.
