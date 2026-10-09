# Native Comparator Manuscript And Figure Integration

Tested integration source was committed/pushed at
`525e27b9b78243793ecd9b3b944514c623c3d126` before actual generation. The
focused/adjacent suite passed 69 tests in 4.03 seconds, including a temporary
real figure and exact parent-body restoration. No new bootstrap, inference,
conversion, scoring, resource admission or benchmark retry occurred.

## Actual Outputs

The [new manuscript](native_qfo_comparator_manuscript_20261009_v1.md) is 82,500
bytes, SHA256
`cb6741dcc53f5b67bb20deccdefe5f165c20f001385b227f697e0f19571d2673`.
Its new Methods and Results insertions integrate the checked retrospective
SwissTrees analysis. Removing precisely those insertions restores the exact
terminal-failure parent manuscript. Earlier evidence, failure descriptions,
endpoint tables, citation IDs, defaults and archives remain unchanged.

The [six-panel comparison figure](native_qfo_comparator_integration_20261009_v1/native_comparator_sensitivity.pdf)
uses percentage points, whereas the
[48-row plotted-endpoint table](native_qfo_comparator_integration_20261009_v1/endpoints.tsv)
and manuscript F1 table retain raw 0-to-1 units. Full OrthoFinder comparisons
are above sequence-only checkpoint comparisons. All eight cells appear in
both rows; four unavailable cells have explicit reasons rather than zero dots.
Thick lines are nominal 95% intervals and thin lines are the separate
48-endpoint adjusted intervals. The earlier 24/42-endpoint analyses keep
their own scopes.

The [generation manifest](native_qfo_comparator_integration_20261009_v1/manifest.json)
is 3,215 bytes, SHA256
`5e79563ab71fac550c358aa20cce4ed03263577602e97e8c18daa503e25a1128`.
Its initial `visual_reviewed`/`manuscript_rendered` flags remain false because
they record generation-time status; later actual checks are in the separate
[execution/readback receipt](native_qfo_comparator_integration_execution_20261009_v1.json),
49,836 bytes, SHA256
`fc0e811c182ead8bfa9eb811a00116271e6de7c1649cf1b4386c5889a481a8e7`.
This receipt retains six actually executed commands, their code/results and
20 output identities; it is not a promise of future checks.

## Numerical And Rendering Checks

Independent readback imports no integration helpers. All 48 plotted-table
rows, 24 estimated endpoints, 24 missing endpoints, eight manuscript F1 rows,
and exact parent-body restoration passed. Direct examination of PDF vector
coordinates checked all 48 interval segments and all 24 point markers against
the source result with the actual axis transforms, colors and line widths.
All 85 figure PDF text spans are inside page bounds. The actual figure PNG
and independently decoded one-page PDF preview were visually inspected;
labels, intervals, unavailable statuses and caption were readable without
observed clipping or incoherent overlap.

The [HTML review](native_qfo_comparator_manuscript_20261009_v1_review.html) and
[asset receipt](native_qfo_comparator_manuscript_20261009_v1_assets.json)
resolve 19 citation IDs against the exact parent-bound bibliography, which
contains 38 entries, not 19. All 121 local link occurrences and 117 distinct
targets were checked. The new figure was correctly recorded as untracked at
render time; later staging does not rewrite that historical receipt.
Pandoc parse/render stderr is empty.

The [actual printed manuscript](native_qfo_comparator_manuscript_20261009_v1_print/document.pdf)
has 24 pages, 440,373 bytes and SHA256
`4ca152292c73301bd3a745f1abc8e2626d2b57e42927462039931f3e1230bb36`.
Its [print receipt](native_qfo_comparator_manuscript_20261009_v1_print/print.json)
records the zero exit and fresh private browser profile. The
[layout check](native_qfo_comparator_manuscript_20261009_v1_layout/report.json)
finds no page-bound violations. All 24 printed F1 numeric cells and eight
family-outcome strings match the source result; selected protocol/scope phrases
are present. Pages 3-5 and 9-12, covering the insertions and adjacent content,
were actually viewed and found readable. The other 17 pages were not manually
viewed in this revision: no full-document visual certification is claimed.
The figure is linked, not embedded in the manuscript PDF.

## Scientific Scope And Next Work

The [scientific result and claim limits](NATIVE_QFO_COMPARATOR_RESULT_20261009.md)
remain unchanged: all four native F1 point estimates are below full
OrthoFinder; reconciled cells show higher precision/lower recall against
the sequence-only checkpoint, with adjusted F1 intervals including zero.
Neither equivalence, a new default, an initial-HMM causal effect nor general
superiority follows. These conditional intervals cover SwissTrees only.

The existing three-cell fixed-bin and model-distance reports describe
profile-off cells. They cannot estimate the newly admitted profile-on contrast.
A separately frozen descriptive extension can project that cell into the
unchanged existing bins without new raw-feature extraction or outcome-based
cutoffs. This is original goal 4.3 work; it does not resolve literal fragment,
ancestral-history or complete architecture truth. Other endpoint uncertainty
and historical provenance limitations stay explicit. The full publication goal
remains active; no new archive candidate or publication-readiness claim follows.
