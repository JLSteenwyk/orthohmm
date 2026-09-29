# OrthoBench Descriptive Figure Evidence Bundle

Built from committed source `72af07f885b465c8c6241ca6fc1fed695f13824b`
using the explicit `ob-complete-strata` audit. The bundle contains nine
files, 1,014,348 payload bytes, one panel and three outputs (PNG/PDF/SVG).
Its manifest SHA-256 is
`133ead009869f6a56121f4eca03ac8d6d1f5e39e658cad08cb7b322d6f120e31`.

A copy outside the checkout passed the bundled standard-library verifier
under `/usr/bin/python3 -I -B`, with working directory `/tmp`. The
[relocation receipt](ob_complete_figure_bundle_relocation_20260928.json)
records the actual external directory, exact command, zero exit code and
matching inventory. No scientific inference or plotting was rerun.

The extended manuscript now includes the supplementary figure, its source
table and bounded caption. The reproduction guide links the export and
figure. Existing rendered manuscript snapshots have not been regenerated
and predate that insertion; the concise main-text draft is unchanged.

## Reproduce

Choose a fresh output directory from the repository root:

```sh
python -B -m benchmark_tools.bundle_publication_figures build --repo . --revision 72af07f885b465c8c6241ca6fc1fed695f13824b --audit benchmark_tools/results/ob_complete_figure_integrity_20260928.json --output /tmp/ob-complete-figure-bundle
python -I -B /tmp/ob-complete-figure-bundle/benchmark_tools/bundle_publication_figures.py verify /tmp/ob-complete-figure-bundle
```

This verifies direct figure evidence and explicit relocation mappings. It
does not include all imported plotting helpers, compiled dependencies,
transitive raw data or native inference. It is not a standalone executable
analysis, redistribution-rights determination or complete publication
archive. The score report retains historical absolute provenance paths;
the verifier treats these as mapped identities rather than opening them.
No new uncertainty or generalization claim is added.
