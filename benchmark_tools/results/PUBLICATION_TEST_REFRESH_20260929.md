# Full Unit Regression And Parser Follow-Up

At `cc3df68a1a946e0a4b3b1bc00b15a079472a2f96`, the full `tests/unit`
suite completed with **11,885 passed, 12 skipped, 22 warnings**, exit zero.
Pytest reported 456.21 seconds; the parent execution window was 467.09
seconds. These are test durations, not inference benchmarks. The previous
full run at `3eb446e` remains historical: 183 source/test files changed
between the two revisions. No fix or retry was needed for this full run.

The [machine-readable receipt](publication_unit_refresh_20260929.json)
records command, interpreter, revision, JUnit totals, log/JUnit hashes and
every skip reason. Tracked status for `tests/unit`, `benchmark_tools/*.py`
and `orthohmm` was clean before and after. Unrelated files were preserved.

Ten skips remain opt-in installed-native tests. Two additional skips were
module-collection skips in the Selectome inspectors because the shared
interpreter lacks `sqlglot`; they do not represent only two skipped test
cases. A separate targeted run exposed the already retained private
`sqlglot` 30.20.0 through `PYTHONPATH` and passed all **13 tests**, with no
skips. Its JUnit report is separately identified in the receipt. No package
was installed and the shared interpreter was not modified. This is not a
single combined full-suite pass in the parser-enabled environment.

The 22 warnings are eleven each for existing invalid escape sequences in
frozen `parser.py:24` and `writer.py:57` exercised by archive verification.
No frozen scientific code was edited to silence them.

## Executed Commands

```bash
python -m pytest tests/unit -q \
  --junitxml=benchmarks/work/publication_unit_refresh_20260929/junit.xml
env PYTHONPATH="$PWD/benchmarks/work/treefam_selectome_search_20260928/parser_lib" \
  python -m pytest -q tests/unit/test_inspect_selectome_tf7a.py \
  tests/unit/test_inspect_selectome_tf7ab.py \
  --junitxml=benchmarks/work/publication_unit_refresh_20260929/selectome.xml
```

Use new report paths for another run. This validates unit-level behavior at
the recorded revision, not independent biological accuracy, controlled timing,
native integration at this revision, rights clearance or archive restoration.
The original scientific failures and remaining publication requirements stay
open.
