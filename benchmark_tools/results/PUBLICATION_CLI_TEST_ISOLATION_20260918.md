# CLI Test Isolation And Full Unit Check

The full unit run at source b3217b0 reported5512passed,9skipped,4failed
in114.40s. All failures were in `test_entrypoint.py`: the shell could not
find an installed `orthohmm` command. This is retained as a failed run,
not silently replaced by the subsequent result.

Inspection also found that its old success-path test ran inference directly
in tracked `tests/samples/`, without checking output content. That path was
not reached in this failed run. The unit tests now invoke the local module
with the running Python interpreter for help, no-argument and missing-input
behavior, without a shell or PATH-dependent wrapper. Successful parser-to-
executor dispatch is tested with temporary FASTAs and a mocked executor.
The existing zero exit status for missing input is preserved, not endorsed
or changed by this test-only work.

A separate integration test executes `python -m orthohmm` end-to-end with
built-in search and Leiden, using temporary copies of the sample FASTAs.
It checks successful completion, nonempty groups, exact input-gene partition
coverage with no duplicated or missing members, and unchanged source/copied
FASTA bytes. It does not assert that this small fixture proves biological
accuracy or general determinism. Installed console-script packaging remains
a separate check; local-module execution is not claimed to cover it.

## Verification

- `python3 -m pytest tests/unit/test_entrypoint.py -q`:4passed in0.57s.
- `python3 -m pytest tests/unit -q --junitxml=benchmarks/work/publication_unit_entrypoint_fixed_20260918.xml`:
  **5516passed,9skipped in117.63s**.
- `python3 -m pytest tests/integration/test_module_cli.py -q --junitxml=benchmarks/work/publication_module_cli_20260918.xml`:
  **1passed in3.92s**. The other integration tests were not rerun in this check.

Retained JUnit files in `benchmarks/work/`, SHA-256:

| File | SHA-256 |
| --- | --- |
| publication_unit_b3217b0_20260918.xml | c95db9128b24c31a0a9421c8d675d440714fc8b3ee032509958f4f590c3c1cd4 |
| publication_unit_entrypoint_fixed_20260918.xml | cd807662eb6868d214b2e71e0d63699c5f396015b1b4d03a9bb563489cc0c72b |
| publication_module_cli_20260918.xml | ca9f9fc5e66e00de253e87e11496a03e8ac993fb4add977ae4327b1971e24d72 |

Only tests changed. No inference code, frozen executor, benchmark score,
default parameter or active timing collector was modified. Existing unrelated
sample-output changes were left intact.
