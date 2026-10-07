# Four Cell Manuscript Render Bibliography Selection

The first v4 HTML render failed before producing either HTML or its asset
receipt. The [retained failure](four_cell_main_render_selection_failure_20261007_v1.json)
records actual exit 1 and the missing citation `iqtree3_2026`. The selected
historical September CSL v5 does not contain that citation. This was a
bibliography-selection error in the handoff, not a renderer or manuscript
defect. The scientific source and all historical evidence remain unchanged.

The [actual parent v3 render receipt](publication_main_render_20261007_v3.json)
binds `publication_bibliography_20261007_v1.csl.json`, 50,517 bytes, SHA256
`a8b74a33a5244a155779567aaf838895f03645c57e9e46be16b30105ce56f392`.
Fresh preflight checks that exact binding, all 19 unchanged manuscript citation
IDs, unique bibliography IDs and absence of missing citations. Every old v5
entry is unchanged; exactly `iqtree3_2026` was added in the already-retained
October bibliography. No metadata download, new citation or source edit is
needed. Keep the failed selection explicit; do not silently relabel it success.

Commit/push this prospective selection correction and failure before one new
render using the unchanged renderer, the exact v4 source and this already
verified October bibliography. Fresh destinations:
`PUBLICATION_MAIN_TEXT_20261007_v4_checked.html` and
`publication_main_review_20261007_v4_checked_assets.json`.
On success use unchanged print/review tools in fresh
`publication_main_review_20261007_v4_checked_print` and
`publication_main_review_20261007_v4_checked_pdf_review` namespaces.
Retain actual commands, outcomes, bindings and full-document visual review.
No automatic loop, new scientific admission or modification of any existing
source, bibliography, successful or failed artifact is authorized by this
correction. Earlier references to rendering v4 with CSLv5 are superseded only
for the bibliography selection and output namespace.
