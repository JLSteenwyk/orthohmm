# OrthoHMM Study Candidate 2026.10.04-rc2

This local candidate adds the main OrthoBench/SwissTrees comparison and
factorial statistical replay routes to rc1's source, manuscript, figures,
resource reporting and native-execution assets. Scientific settings, scores,
intervals and the32-page review remain unchanged. Software remains0.5.0.
No public release, DOI, submission or publication readiness is certified.

## Start Here

Retain the archive and `PACKAGE_INDEX.json` digests outside the package. Use
the copied standard-library reader, an unused output path and external anchors:

```sh
python3 -I -S -B /trusted/bundle_publication_package.py restore /candidate.tar.gz \
  --archive-sha256 ARCHIVE_SHA256 --manifest-sha256 PACKAGE_INDEX_SHA256 \
  --output /fresh/package
python3 -I -S -B /fresh/package/bundle_publication_package.py verify /fresh/package \
  --manifest-sha256 PACKAGE_INDEX_SHA256
```

The reader verifies stream bytes before extraction and exact inventories
afterward; it never starts scientific jobs. Keep all generated replay output
outside the immutable package. No Git, checkout imports, scientific packages
or workstation data paths are required for this integrity verification.

| Location | Current scope |
| --- | --- |
| `review/document.pdf`, `figure-replay/` | Unchanged32-page review and original presentation inputs |
| `archives/`, `anchors/` | Exact rc1 handoff, final resource reporting and native assets/indexes |
| `statistics/` | Main comparison/factorial count files, workers, imported numerical helpers, protocols and NumPy requirement |
| `evidence/` | Prior scoped evidence plus full-goal audit and current all-tool/Three Kingdoms provenance |
| `evidence/PUBLICATION_PACKAGE_20261004.md` | Preserved rc1 outer guide for unchanged nested components |

## Main Statistical Replay

Use a separate Python environment with NumPy2.2.6. The executed original
arithmetic routes used Python3.10.13; the current packaging readback records
its actual interpreter. NumPy is not bundled. This environment is statistical,
not the scientific/native runtime or an OS-hermetic installation.

```sh
python3 -m venv /fresh/statistics-env
/fresh/statistics-env/bin/python -m pip install -r /fresh/package/statistics/requirements.txt
```

Set `PACKAGE` to the restored root, `ANALYSIS_PYTHON` to that environment's
Python, and `FRESH` to an existing writable directory outside the package.
All four output paths below must be absent. One numerical thread is used to
avoid unnecessary contention; these executions are not timing benchmarks.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$ANALYSIS_PYTHON" -I -B "$PACKAGE/statistics/benchmark_tools/check_ob_complete_uncertainty.py" \
  --portable "$PACKAGE/statistics/benchmark_tools/results/ob_complete_portable_statistics_20260928.json" \
  --sha256 cdb30ff987100843f3127a39a75eb3649123a136229b8b41e6ed156209aff30e \
  --output "$FRESH/ob-comparison"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$ANALYSIS_PYTHON" -I -B "$PACKAGE/statistics/benchmark_tools/reproduce_corrected_swiss_comparison.py" \
  --results "$PACKAGE/statistics/benchmark_tools/results/qfo_recovered_swiss_uncertainty_22178.json" \
  --results-sha256 a599749d66433ec211ed1ab0e3a6a2a75eb99cb0dc57abc47c2d267930cbbe34 \
  --retained-counts-only --output "$FRESH/swiss-comparison.json"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$ANALYSIS_PYTHON" -I -B "$PACKAGE/statistics/benchmark_tools/reproduce_orthobench_factorial.py" \
  --worker "$PACKAGE/statistics" --output "$FRESH/ob-factorial"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$ANALYSIS_PYTHON" -I -B "$PACKAGE/statistics/benchmark_tools/reproduce_corrected_factorial.py" \
  --worker "$PACKAGE/statistics" --output "$FRESH/qfo-factorial"
```

The workers recompute the full benchmark statistic from retained family
sufficient statistics and compare against frozen expected results. They
preserve70 OrthoBench families/eight methods/21 endpoints,18 SwissTrees
families/eight methods/24 endpoints, and both eight-cell factorials. The
OrthoBench factorial uses20,000 paired draws; the other three use100,000.
Seeds, multiplicity inventories and tolerances are unchanged. Original
raw-source/reference admission and exchangeability assumptions are not proved
by arithmetic reproduction. No native inference, pair conversion or raw
scoring is executed. Uncertainty for other QfO endpoints remains unresolved.

## Other Reproduction Routes

The preserved rc1 guide supplies presentation replay, final resource-table/
PNG/prose replay, trusted handoff/native-asset restoration, raw OrthoBench
acquisition, private-base preparation and the integrated native controller.
Those unchanged routes retain their original exact anchors and execution
receipts. Do not rerun admitted job22377, completed scientific diagnostics or
any of the27 timing identities solely because this candidate exists.
YGOB and fixed/variable simulation arithmetic remain in the nested handoff's
`arithmetic/` directory with their documented requirements and limitations.

QfO/OrthoBench are primary and development-exposed; Three Kingdoms is a
supplementary BUSCO-restricted endpoint. The all-tool provenance register
distinguishes native pairs, group-derived pairs and checkpoints, inherited
identities and missing historical evidence. Current timing uses the shared
Threadripper, with unknown potentially method-dependent contention. No DGX,
quiet window, unrelated-job disruption or causal efficiency ranking is needed.

## Remaining Boundaries

This adds a direct main-statistics delivery route, not full transitive raw
inference for every competitor. Raw inputs/references, base/bootstrap/OS,
some third-party identity/redistribution evidence and historical per-ablation
costs remain separate. Original TreeFam-A7 family files are unavailable in the
inspected public sources. Several QfO uncertainty and prespecified error-stratum
requirements remain unresolved. Read `evidence/PUBLICATION_REQUIREMENT_AUDIT_20261004.md`
as a dated full-scope assessment, alongside the newer provenance supplements.
Public deposition and submission remain unexecuted. Preserve rc1 unchanged;
this candidate does not retroactively certify it or the full publication goal.
