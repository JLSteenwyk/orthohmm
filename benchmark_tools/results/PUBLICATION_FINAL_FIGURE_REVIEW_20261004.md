# Current Main Text With Complete Figure Appendix

The [new working review](publication_main_with_figures_20261004/document.pdf)
contains32 pages: twelve main-text pages, four caption-guide pages and sixteen
original vector figures. The historical nine-page/no-resource presentation
and its accepted replay remain unchanged. This is a presentation checkpoint,
not a submission-ready or complete executable study release.

The [relative-path, hash-pinned selection](publication_figure_selection_20261004.json)
reuses the fourteen accepted scientific figure captions and adds the existing
all-method OrthoBench descriptive strata and final Threadripper resource
figure. Both additions match their producer manifests. Error-strata captions
retain development exposure, small-bin and descriptor limitations. Resource
captions retain27 reviewed attempts,25 measurements,24 eligible observations,
six complete cells and three unavailable cell summaries; unknown contention
does not justify an isolated efficiency ranking.

The [assembly manifest](publication_main_with_figures_20261004/assembly.json)
records exact geometry, text and pixels for all28 original source pages.
Seven main-text figure links now target embedded figures, all sixteen guide
entries have internal destinations, and eighteen bookmarks are present.
Fifty-eight original nonfigure file-URI actions retain their historical
locations; those links do not become portable data access.

Actually inspect the four new guide pages and two added figure pages: no
observed clipping or overlap. The [review/test receipt](publication_final_figure_review_20261004.json)
distinguishes that manual scope from exact preservation of the prior main
and fourteen figure pages. No fresh review of every original page is claimed.
The assembly's original visual-review flag stays false rather than being
silently relabeled after inspection.

Thirteen focused cases pass in6.46s. They exercise normal and optimized
isolated CLI replay, link destinations/actions, preserved source pages,
existing-output refusal, changed bytes, unsafe paths, duplicate/ambiguous
inventory, wrong page counts, overflowing captions and symlink escapes.
The actual sixteen-figure inventory also replays from a copied, input-only
tree with isolated optimized Python, without checkout imports or Git;
all32 pages match canonical geometry, text, pixels and link actions.
The copied fixture is under ignored benchmark work, not a delivered archive.

Rebuild from the repository root into a new directory:

```sh
python -B -m benchmark_tools.assemble_publication_review --root . --selection benchmark_tools/results/publication_figure_selection_20261004.json --selection-sha256 4f21c1ea895d1d2fa26e345afbeb18aa961a3df6d37a45b8985474ce33de8f0e --output benchmark_tools/results/publication_figure_review_replay
```

This requires the pinned PDF inputs and PyMuPDF1.27.2.3; it does not rerun
inference or scoring. The inventory builder is separately retained as
`prepare_final_figure_selection_20261004.py`. An eventual versioned archive
must include this selected presentation alongside the current source,
main-text review, final reporting and retained native reproduction assets.
