# Complete Profile-Stratum Manuscript Integration

The integration source was committed/pushed at `0be0d844` after 156 tests
passed. Actual generation and independent numerical readback exited zero.
The [new manuscript](native_qfo_four_cell_strata_manuscript_20261009_v1.md)
contains every one of the 23 fixed profile-bin effects, including five empty
bins, and the two changed-family count records. Removing just its new Methods
and Results sections restores the exact comparator manuscript. Older three-cell
figures/descriptions remain explicitly scoped to their original cells.

## Scientific Checks

The manuscript is 87,595 bytes, SHA256
`60af8bf7ead43a3860eee6c923fc793895a8d8ae8e5564a30b0fe0daef131efc`.
The [generation manifest](native_qfo_four_cell_strata_integration_20261009_v1/manifest.json)
is 2,861 bytes, SHA256
`652f0a593ec8966d55158674a94c47c01a40ad53eaa835a865313d5ffd44ed7a`.
Its generation-time render/visual flags remain false; the later actual reviews
are recorded separately, not retrospectively rewritten.

Independent readback imports no integration helpers and verifies all 23
[profile TSV rows](native_qfo_four_cell_strata_integration_20261009_v1/profile_bins.tsv),
all 69 raw-to-human metric cells, NA status, exact parent restoration,
CASP/GH14 integer count changes and direct input/output/source identities.
The [scientific result and scope](NATIVE_QFO_FOUR_CELL_STRATA_RESULT_20261009.md)
are unchanged: overall F1 decreases by 0.320215 percentage points; small
positive lower-entropy/higher-model-distance bins reflect CASP false-positive
removals, while complementary negative bins contain the GH14 true-positive
loss. These are overlapping descriptive associations, not subgroup confidence
intervals, a total-HMM effect, a causal edge or a new default.

## Actual Rendering Recovery

The first HTML/PDF print succeeded technically and passed page-bound checks,
but visual review of pages 13-14 found the recall column clipped. The
[initial HTML](native_qfo_four_cell_strata_manuscript_20261009_v1_review.html),
[initial print](native_qfo_four_cell_strata_manuscript_20261009_v1_print/document.pdf)
and their receipts remain unchanged and are not labeled visually complete.

The separately tested scoped rendering correction was committed/pushed at
`6b393fee` after 166 focused/adjacent tests passed. It verifies the exact target
HTML structure and adds print CSS to that one table only. Removing the inserted
style restores every byte of the original scientific HTML. It changes neither
the Markdown manuscript nor any score, table value, citation, link or old global
print header. The new [HTML review](native_qfo_four_cell_strata_manuscript_20261009_v1_review_v2.html)
and [asset receipt](native_qfo_four_cell_strata_manuscript_20261009_v1_assets_v2.json)
preserve the original executed Pandoc command as inherited provenance, not a
claim that Pandoc ran again. The original render resolved 19 citation IDs,
128 local link occurrences and 123 targets with empty stderr.

The [corrected print](native_qfo_four_cell_strata_manuscript_20261009_v1_print_v2/document.pdf)
has 25 pages, 507,887 bytes, SHA256
`f267d80682aa25de8e31bb65fc59c0b58970cd8382ca18e5da4d0ce8bd9d8d59`.
An independent coordinate readback checks every printed suite/bin label,
family size and all 69 metric-or-NA cells against the scientific source.
All numeric cells fit their actual columns and page bounds; all six columns
are visible. A first inline readback stopped on an HTMLParser helper-name
collision before PDF checking; its code/error is retained, and a separately
corrected inline reader passed without changing any input or PDF.

The [layout review](native_qfo_four_cell_strata_manuscript_20261009_v1_layout_v2/report.json)
finds no page-bound violations. Actual pages 4-6 and 12-15 were viewed and
found readable without observed clipping or incoherent overlap. The other
18 pages were not manually viewed; this is not full-document visual
certification. The [combined execution record](native_qfo_four_cell_strata_integration_execution_20261009_v1.json)
retains ten actual command/code/tool outcomes, both test runs, all 29 selected
artifact identities, the initial visual failure and precise later review scope.
Private browser profiles are not publication artifacts or staged Git content.

## Remaining Scope

The report and integration resolve the missing profile-on fixed-bin projection
and manuscript coverage. They do not resolve literal fragment/complete-domain
truth, known ancestral histories, missing factorial scores or paired uncertainty
for every QfO endpoint. No inference, benchmark scoring, feature extraction,
bootstrap, accuracy/resource admission or archive generation occurred here.
Old source/evidence/drafts/archives remain unchanged. The full seven-part goal
remains active; publication readiness and general superiority over full
OrthoFinder remain unproved.
