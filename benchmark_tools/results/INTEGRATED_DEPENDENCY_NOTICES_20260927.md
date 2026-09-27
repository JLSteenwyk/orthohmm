# Integrated Workflow Dependency Notices

Collected exact embedded notice candidates from the wheel artifacts selected
by the two completed installation stages of job 22337. The native job is still
running; this work does not admit its scientific results or change its files.
The [inference inventory](integrated_inference_notices_20260927.json) and
[reader inventory](integrated_reader_notices_20260927.json) are recomputed from
hash-pinned pip reports and wheel metadata before export.

| Distribution | Version | Notice candidates | Native members |
| --- | --- | ---: | ---: |
| Biopython | 1.87 | 1 | 13 |
| DendroPy | 5.1.0 | 3 | 0 |
| igraph | 1.0.0 | 1 | 4 |
| leidenalg | 0.11.0 | 1 | 3 |
| llvmlite | 0.49.0 | 2 | 1 |
| numba | 0.67.0 | 2 | 14 |
| NumPy | 2.2.6 | 4 | 22 |
| OrthoHMM | 0.5.0 | 1 | 3 |
| pip | 26.2.1 | 44 | 0 |
| python-igraph | 1.0.0 | 1 | 0 |
| setuptools | 83.0.0 | 19 | 0 |
| texttable | 1.7.0 | 1 | 0 |

The environments contain 11 and five distributions, respectively. Shared
wheel digests are counted once: 12 distinct artifacts, 80 notice candidates
(488,560 bytes) and 60 native-library filename candidates. All declared
license-file fields resolve uniquely in these artifacts. Provider declarations
remain verbatim in the inventories; they are not compatibility determinations.
Unlike the older development wheelhouse inventory, this set contains the
actual recovery Leiden 0.11.0 wheel, frozen scientific OrthoHMM wheel and
patched reader Biopython 1.87 wheel.

## Export And Verification

The [export receipt](integrated_notice_export_20260927.json) records a local
notice-only directory, `benchmarks/work/publication_integrated_notices_20260927`.
Files are separated by full wheel digest and retain their archive member paths.
Each member's bytes are checked before writing and after export. The exporter
rejects unsafe paths, nonregular members, conflicting duplicate inventories,
omitted/extra candidates and changed payloads. No wheel code or binaries are
exported or loaded.

The 73,506-byte `NOTICE_INDEX.json` has SHA256
`e14681f6d8b4dba59bd4d198df4d7e9ceb41f9dbe3dac62cff014b986de1d6f5`.
All 80 texts and the identical index verify after copying to a different root;
see [relocation receipt](integrated_notice_relocation_20260927.json).
Twenty-eight focused inventory/export tests pass.

```bash
python -m benchmark_tools.export_dependency_notices \
  --inventory benchmark_tools/results/integrated_inference_notices_20260927.json \
  --inventory benchmark_tools/results/integrated_reader_notices_20260927.json \
  --output /fresh/integrated-notices --receipt /fresh/integrated-notices.json
```

This collects review material; it does not close the redistribution review.
Filename heuristics can miss notices, and AUTHORS files need not be licenses.
Native-member filenames do not reveal every statically linked component or
establish a complete component-to-notice mapping. External tools, OS libraries,
containers and datasets remain outside scope. No legal compatibility verdict,
public deposition, runtime redistribution or publication-readiness claim is made.
