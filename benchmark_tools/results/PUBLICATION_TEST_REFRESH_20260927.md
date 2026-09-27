# Current Unit Regression Result

At source revision `9914060112f57b79e75401aeb907c3a6d6b2b05f`, the full
`tests/unit` suite completed with exit zero:

- 10,729 passed.
- 10 skipped, all opt-in installed native-runtime checks.
- 22 warnings: eleven repetitions each of existing invalid-escape warnings
  in frozen `parser.py:24` and `writer.py:57` during source-archive tests.
- Terminal elapsed time: 315.96 seconds; this is not a benchmark runtime.

The [machine-readable receipt](publication_unit_refresh_20260927.json)
records the tested revision, exact command/interpreter, pytest version,
JUnit checksum, suite metadata and all ten skipped test names/reasons.
The full JUnit report is retained at
`benchmarks/work/publication_unit_refresh_20260927.xml`.

Tracked source status was clean for `tests/unit`, `benchmark_tools/*.py`
and `orthohmm` before and after the run. Unrelated sample changes remain
untouched. No code fix or retry was needed.

The skipped checks concern installed legacy BLAST, OrthoMCL/BioPerl,
database normalization, checkpoint admission and native staging. They must
not be counted as current native-integration coverage. This pass does not
establish held-out biological accuracy, dedicated timing, release portability,
third-party rights or publication readiness. Historical scientific sources
remain unchanged despite their syntax warnings.
