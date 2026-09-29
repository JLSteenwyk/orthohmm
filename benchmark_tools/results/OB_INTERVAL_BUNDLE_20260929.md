# Relocated OrthoBench Interval Figure

The new explicit `ob-complete-intervals` audit scope checks the complete
21-endpoint design, units, pinned input and retained figure bytes. Forty
focused audit, plotting and export tests pass. The audit was committed and
pushed as `18b9f991` before export from that exact revision.

The existing bundle exporter produced 11 files (357,535 bytes), one panel and
four outputs, plus its relocation manifest. A copy under
`/tmp/orthohmm-ob-interval-relocation-20260929/bundle` verified using the copied
standard-library verifier with `/usr/bin/python3 -I -B` from `/tmp`.

More than file verification was tested: the copied plotter was executed with
the copied result and helper using Python isolated mode from `/tmp`. Both
project modules were explicitly checked to resolve inside the copied bundle.
Regenerated `endpoints.tsv` and PNG bytes match the retained versions exactly.
PDF and SVG were generated but are not claimed byte-identical because their
metadata can vary. The original visual review remains applicable to the
byte-identical PNG. No new statistical calculation was performed.

- [Figure integrity audit](ob_interval_figure_integrity_20260929.json)
- [Relocation and actual render commands](ob_interval_bundle_relocation_20260929.json)
- [Versioned archive receipt](ob_interval_bundle_archive_20260929.json)

The local archive is
`benchmarks/work/ob_interval_bundle_18b9f991.tar.gz`. All 12 regular-file members
were reopened and compared byte-for-byte with the exported bundle. No archive
member was extracted or executed during that archive check.

```sh
python -B -m benchmark_tools.bundle_publication_figures build \
  --repo . --revision 18b9f991 \
  --audit benchmark_tools/results/ob_interval_figure_integrity_20260929.json \
  --output /tmp/ob-interval-bundle-new
python -I -B /tmp/ob-interval-bundle-new/benchmark_tools/bundle_publication_figures.py \
  verify /tmp/ob-interval-bundle-new
```

This is same-host relocated plotting with existing scientific Python packages,
not a hermetic environment reconstruction, native inference reproduction,
raw-data archive, third-party rights clearance or completed publication release.
An initial shell-quoting error prevented a preliminary relocation command from
executing; it created no files. The retained successful command is in the
relocation receipt. No benchmark rerun, scheduler change or timing admission.
