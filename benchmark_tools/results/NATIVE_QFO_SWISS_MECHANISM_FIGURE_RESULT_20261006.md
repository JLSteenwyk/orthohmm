# Native SwissTrees Mechanism Companion Figure

## Delivered Artifact

The final four-panel [PNG](native_qfo_swiss_mechanism_figure_20261006_v3/native_swiss_mechanism.png),
[PDF](native_qfo_swiss_mechanism_figure_20261006_v3/native_swiss_mechanism.pdf)
and [SVG](native_qfo_swiss_mechanism_figure_20261006_v3/native_swiss_mechanism.svg)
join the already independently checked native P0/C0/R0 and P0/C0/R1
diagnostics. Initial HMM search is on in both cells. These are native ablation
cells, not selected-default/high-sensitivity competitor results. The figure
is a manuscript companion; it is not yet in the historical manuscript or
archive and does not establish publication readiness.

Panel A includes all 2,023 changed SwissTrees pairs, not chosen examples:
334 removed TP have 60 no-hit, 6 one-direction and 268 bidirectional cases;
1,689 removed FP have 479, 20 and 1,190 cases, respectively. Selected directed
score multisets are identical between the original search checkpoints. This
does not imply whole-hit-set equality or orthology confidence. Absence of a
significant direct hit does not distinguish prefilter loss from scoring or
indirect graph recruitment.

Panel B reports the same complete cohort: 280 removed TP remain within a
final RootHOG and 54 cross final roots; all 1,689 removed FP remain within a
final RootHOG. This separates within-group pair exclusion from group splitting.
The previously checked reconciliation annotations locate those exclusions;
neither annotations nor this figure prove biological event correctness.

Panels C and D show native F1 and macro precision/recall differences in all
18 families, frozen lower/higher-entropy bins of nine each, seven short-relative
families and the other 11. F1 is the harmonic mean of macro precision and
recall, not pooled-pair F1 or mean family F1. Bins are input-only, unchanged,
overlapping and development-exposed. There are no subgroup intervals,
significance tests or new uncertainty draws. Entropy is not divergence or
local low complexity; relative shortness is not fragmentation. The full-family
F1 shift is +10.039 percentage points, with precision +30.521 and recall
-6.533 points; this figure does not replace the previously admitted conditional
SwissTrees interval, which includes zero for the F1 difference.

All values retain full machine precision in three generated ten-row tables:
[stages](native_qfo_swiss_mechanism_figure_20261006_v3/stages.tsv),
[strata](native_qfo_swiss_mechanism_figure_20261006_v3/strata.tsv) and
[differences](native_qfo_swiss_mechanism_figure_20261006_v3/differences.tsv).
The three versions have byte-identical corresponding tables, verified by
direct byte comparison and a retained regression test. Only layout changes.

## Evidence And Checks

The renderer binds the actual search report/integer-key readback,
reconciliation report/saved-Newick readback and sequence-strata report/stdlib
readback. Current sources, report pins, readback scope, exact case identities,
RootHOG flags, common native admissions, complete counts and plotted arithmetic
are checked. Direct bound files are checked before and after rendering.
This reuses scientific evidence; no 181-million-row array scan, raw-scoring
replay, inference, tree reconstruction, reconciliation or transitive accuracy
admission is repeated simply to make a figure.

The [independent final readback](native_qfo_swiss_mechanism_figure_readback_20261006_v3.json)
imports no renderer or aggregation helper. It derives and compares all 30 TSV
rows and exact headers against the original bound reports, checks the six
output identities, SVG scope/data labels and four nonblank raster panel crops
in both the PNG and PDF-derived bitmap. PNG dimensions are 2,640 by 1,760.
Pixel-color/crop checks establish visible assets, not mathematical validation
of plot geometry or biology.

Both the actual final PNG and [PDF bitmap](native_qfo_swiss_mechanism_pdf_preview_20261006_v3.png)
were visually inspected with `view_image`; labels, titles, legends, axes,
points, bars and limitations are visible and do not overlap. This is agent
visual inspection, not human peer review. The automatic receipt deliberately
does not set `human_visual_review_complete`. A separate actual-Matplotlib test
draws the figure and checks every legend's display-coordinate bounding box
against its plot area and title, as well as patch counts and all three formats.
Its success is distinct from the independent raster checks.

Retain the initial 28-test success and failed joined 40-case receipt:
39 passed and the actual geometry test failed. Version 1 has visibly overlapping
lower-panel legends. Version 2 improved their placement, but strict bounding
boxes still overlap every plot area slightly. Version 3 anchors each legend
outside its axes and clears its title, without relaxing that test. Historical
versions and manifests remain unchanged; their older renderer source hashes
are not treated as current-source verification.

Retain the actual [failed version-2 reader record](native_qfo_swiss_mechanism_readback_failure_20261006.json).
A local RootHOG `Counter` shadowed the directory argument, causing exit 1
before the requested JSON existed. Rename the local variable only; the fixed
version-2 reader succeeds under a new output path. Its original pre-fix source
hash is recorded, but those source bytes were not separately archived.
The final version-3 readback is a separate actual invocation, not a renamed
old result. No native inference or scoring was retried for this repair.

The final joined regression suite passes all 394 cases in 17.81s, without
failures, errors or skips. The new renderer/reader contribute 42 cases; the
suite also covers existing
sequence-strata, search support, reconciliation, transitions, ordinary/recovered
counts, Swiss uncertainty binding and scientific score export. Its final XML
is `native_qfo_swiss_mechanism_complete_joined_tests_20261006_v3_final.xml`.
The initial 392-case success is also retained; two final artifact-history and
actual-receipt binding checks extend the suite. Tests do not establish
whole-study completion or supply missing native cells.

Rendering uses sanitized Python 3.12.3, Matplotlib 3.10.8 and single-threaded
numeric-library settings. The final GNU-time receipt records 1.43s,
112,728 KiB maximum RSS, zero swaps and exit 0. These are shared-host
postprocessing observations, not inference timing or an efficiency comparison.
CPU, memory-bandwidth and I/O contention have unknown, potentially tool-dependent
effects. Capacity at this launch: roughly 626 GiB available RAM and 9.5 TiB
disk; host swap remains almost full. No unrelated workload was modified.
The failed R1 native timing remains ineligible and no timing is plotted.

Final renderer source SHA256:
`26b1fd50424d47914208417e1cb2a2a3a9c4d0222572270f24ad89d80b78c180`.
Independent reader source SHA256:
`84202b96a761b93aac8c9dbdebd9faedf969f5a9183a9f6c4b3337e1f7c6417d`.
Final manifest: 9,002 bytes, SHA256
`c3e96b8a2aecc8d9926c19875e3a321c4edb1173d522b012301152da137876fe`.
Final readback: 10,213 bytes, SHA256
`51dc6135ddc367c1dff32361f70e22d365c667c417274d0d9f545b35f2d98625`.
The manifest records exact input/output/evidence identities and narrow scope.

## Reproduction

From the repository root, with the bound local diagnostic artifacts available:

```bash
env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/release_alert_refresh_20261001/venv/bin/python -B \
  benchmark_tools/plot_native_qfo_swiss_mechanism.py \
  --search-report benchmark_tools/results/native_qfo_swiss_search_support_20261006_v1.json \
  787bde01ffacd60547fb290bdbb1d04242230ed05ce77d25d270227310252d5a \
  --search-readback benchmark_tools/results/native_qfo_swiss_search_support_code_readback_20261006.json \
  8ffc58fe192fa060de4e068161e72b8ab626066446608e5c0011925e62c67ede \
  --reconciliation-report benchmark_tools/results/native_qfo_swiss_reconciliation_trace_20261006_v1.json \
  53c692cf7dab9affac7b0240cf5191d2530bd4369b66f4165e8acb811391dacb \
  --reconciliation-readback benchmark_tools/results/native_qfo_swiss_reconciliation_newick_readback_20261006.json \
  ee83382000bb48b95323209a197e7edb44de918bad6c6b8a1fbf8a1ea6e38835 \
  --strata-report benchmark_tools/results/native_qfo_swiss_sequence_strata_20261006_v1/report.json \
  4650f614d2fdeebd65cbcbf4999c61c83cc549560d2738bfbb5502f686f660e3 \
  --strata-readback benchmark_tools/results/native_qfo_swiss_sequence_strata_readback_20261006.json \
  8d0356bc0e159d204939271578e23a6a5b6efcf6198b75a18c0f3d38e0e67c0a \
  --output benchmark_tools/results/native_qfo_swiss_mechanism_replay

pdftoppm -png -singlefile -scale-to 2200 \
  benchmark_tools/results/native_qfo_swiss_mechanism_replay/native_swiss_mechanism.pdf \
  benchmark_tools/results/native_qfo_swiss_mechanism_replay_preview

env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/release_alert_refresh_20261001/venv/bin/python -B \
  benchmark_tools/review_native_qfo_swiss_mechanism.py \
  --root benchmark_tools/results/native_qfo_swiss_mechanism_replay \
  --pdf-preview benchmark_tools/results/native_qfo_swiss_mechanism_replay_preview.png \
  --output benchmark_tools/results/native_qfo_swiss_mechanism_replay_readback.json
```

Use fresh paths; occupied outputs refuse. Existing manifests have local absolute
paths and depend on their actual bound diagnostic sources. This is not an
outside-checkout whole-study reproduction or a portable completed archive.
Use the [prospective figure protocol](NATIVE_QFO_SWISS_MECHANISM_FIGURE_PROTOCOL_20261006.md).

## Remaining Scope

Original job 22444 remained RUNNING at 7:16:30 during this work; original
22445/22450/22451/22452 were dependency-pending. No unfinished native output
was read, original job restarted or successor released. New admissions and
additional score/uncertainty exports must await the actual original terminal
and downstream production gates.

The original TreeFam-A release-7 trees and `treefam2reference.txt` remain
unrecovered. The preceding public-source follow-up found no verified original
archive or mapping; [Sanger's historical page](https://www.sanger.ac.uk/tool/treefam/)
reports that its resource is no longer available there. Historical
[TreeFam Perl tree-format documentation](https://treesoft.sourceforge.net/tf-perl-api7/Tree.html)
describes formats, not the required original data. This is not proof of global
unavailability. No one was contacted; newer releases or documentation were
not substituted for original family labels.

Remaining native cells and matched-search interactions, graph/tree/other error
strata, wider valid uncertainty, independent generalization, provenance,
whole-study reproducibility, manuscript assembly and release/deposition
requirements remain open. Preserve historical manuscript/PDF/archive bytes
and failed timing evidence. No scientific implementation, defaults, endpoints,
accuracy admission or unrelated jobs change. The full publication goal stays
active, with completion unproven; competing analyses are not a quiet-window
or dedicated-host blocker.
