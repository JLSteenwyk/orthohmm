# Reporting Scope Manuscript Review

The October 12 [manuscript](publication_reporting_scope_manuscript_20261012_v1.md)
integrates the approved VGNC/GO/EC, TreeFam and FAS/secondary-summary reporting
amendments. Its source and [claim checklist](PUBLICATION_REPORTING_SCOPE_CLAIMS_20261012.md)
were committed and pushed at `4e7272bb` before this render. Observed scores,
scientific defaults and historical evidence are unchanged. This is reporting
integration, not a new statistical admission or whole-goal completion audit.

## Executed Checks

The existing isolated Python 3.12 review environment executed the six focused
test modules for manuscript rendering, controlled-fragment and profile-table
styling, PDF printing/review and GO/EC arithmetic auditing: **56 passed in
3.44 seconds**, exit zero. The original October 10 manuscript and checklist
retain their checked SHA256 digests, respectively
`cb5076f91b267b9250e7508a51d56bd25af8a015ba0d4db85c13b14443556db1`
and `f9760a0b36f0b59fd620bb29311c86d5cc2d8a7f9a7f4131ecb7c8a8e57d6120`.

Actual execution used the unchanged `render_manuscript_review.render`,
`render_controlled_fragment_review.run`, `print_manuscript_review.print_review`
and `review_manuscript_pdf.review` helpers with separately named outputs.
The environment unset PYTHONPATH/PYTHONHOME/PYTHONUSERBASE and loader overrides,
disabled user-site and bytecode, and set numerical-library threads to one.
The browser sandbox remained enabled. Both commands exited zero.

The [HTML](publication_reporting_scope_manuscript_20261012_v1_review_v2.html)
and [asset record](publication_reporting_scope_manuscript_20261012_v1_assets_v2.json)
resolve 19 citation IDs, 163 local occurrences and 156 unique tracked targets,
with no untracked target. Only the existing profile-table CSS was applied
after Pandoc rendering. This verifies direct inputs, not external URLs,
transitive provenance, redistribution rights or scientific correctness.

The [37-page PDF](publication_reporting_scope_manuscript_20261012_v1_print/document.pdf)
is 914350 bytes, SHA256
`1cce17858ddf1846ecf20ab9424495e21878dfc388479d84dfa648a1166aa458`.
The [print record](publication_reporting_scope_manuscript_20261012_v1_print/print.json)
records the actual browser command. The
[layout report](publication_reporting_scope_manuscript_20261012_v1_layout/report.json)
finds zero text/image page-bounds violations across all 37 pages.

## Actual Visual Scope

The 13 generated page images were actually inspected: pages 13-18, 27-29 and
32-35. These cover both FAS sampling tables, the comparator/terminal tables,
the complete fixed-profile-bin table, shared-host timing limitations,
Discussion, FAS figure/caption, and all eight rows of both new native-bar
tables. Text, labels and table values are readable without visible clipping
or incoherent overlap. Long tables span pages and some identifiers wrap;
journal layout remains a separate editorial step. The FAS figure distinguishes
offset observed scores from the conditional expected-statistic ranges.
VGNC/GO/EC tables label half-widths, denominators and scored-pair coverage.
This is selected-page inspection, not an all-page visual certification.

## Scientific Limits

All 24 native aggregate hashes and values were checked against retained
admissions before drafting the tables. Exact scorer versions establish VGNC
binomial half-widths and GO/EC Student-t half-widths; a field named `stderr`
is not automatically one SEM. Conventional bars do not correct shared-family
dependence or establish paired significance. TreeFam remains descriptive
supporting evidence, FAS ranges conditional, and the six-metric mean a
descriptive secondary summary. Their documented stronger-uncertainty gaps
are nonblocking under the amended goal, not solved science. No benchmark,
inference, timing, search or uncertainty producer was rerun. Earlier reviews
and local rc6 remain unchanged and do not contain this later revision.
