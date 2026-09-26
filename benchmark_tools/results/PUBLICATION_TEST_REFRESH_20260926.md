# Publication Regression Refresh

## Scope and Initial Failure

Started from `359c9ff8a93c81c36a69edc8de195cbbf99f92d5` on the shared
local host, without DGX access. The first full unit invocation returned
9,756 passed, 10 skipped and one failed in 309.77 seconds. Its failure was
`test_method_executor_identity[fastoma]`: the test expected the superseded
`publication_qfo_corrected_fastoma_assessment_v1` directory.

Production commit `926eb714b36c7582aa1d941c5a0ba29002490e7a` had deliberately
selected the replacement `publication_qfo_corrected_fastoma_assessment_v2`
at `0cc0a96c44e377f4e87a1e579ce12012d3478e4e`. The progress ledger records
that replacement for scoring job 22057, and a fresh `git rev-parse HEAD`
in the replacement checkout confirms the same commit. Updated only the
test's directory expectation and added an explicit expected-commit assertion.
Production pins, frozen executors, scientific defaults and scores are unchanged.

The failed JUnit report is retained, not overwritten:
`benchmarks/work/publication_full_unit_20260926.xml`, SHA-256
`405ce5806550f953734c5bdd172f944c9ffefc14e01aa35747ee54b615310016`.
The changed test file SHA-256 is
`1f99019bfec469625926466c9dabbe030bbbd184f933481ab8f5f75320b049c7`.

## Commands

The initial command used the first JUnit filename above. The verification
command differs only in its output filename. Neither full unit command
explicitly sets numerical thread limits.

```sh
/home/bizon/anaconda3/bin/python -m pytest -q tests/unit \
  --junitxml=benchmarks/work/publication_full_unit_verified_20260926.xml
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ORTHOHMM_LEGACY_BLAST_SMOKE=1 /home/bizon/anaconda3/bin/python -m pytest -q \
  tests/integration/test_module_cli.py \
  tests/unit/test_admit_qfo_corrected_bpo.py::test_independent_complete_recheck \
  tests/unit/test_audit_orthomcl_database.py::test_installed_formatdb_normalization_is_detected \
  tests/unit/test_describe_legacy_residue_deletions.py::test_native_query_and_database_delete_o \
  tests/unit/test_orthomcl_bpo_indexes.py::test_installed_native_bpo_parity \
  tests/unit/test_orthomcl_perl_launcher.py::test_installed_perl_probe_has_no_relative_lookup \
  tests/unit/test_prepare_orthomcl_bpo_checkpoint.py::test_installed_checkpoint \
  tests/unit/test_run_qfo_corrected_blast.py::test_installed_legacy_engine_smoke \
  tests/unit/test_stage_orthomcl_native_inputs.py::test_installed_staged_native_inference \
  --junitxml=benchmarks/work/publication_native_checks_20260926.xml
```

## Full Unit Verification

The rerun passed **9,757 tests**, with **10 skips**, in 304.51 seconds.
Parsed JUnit confirms 9,767 total cases, zero failures and zero errors;
its suite duration is 304.132 seconds. SHA-256:
`ef31d8a6b0020b909361db20240a8c526e807c86bb81ce2cde89845fcef6c098`.
The 22 warnings concern invalid escape sequences in the frozen archive's
parser/writer files. No frozen source was edited to silence them.

## Native Checks

All 11 native checks passed in 23.94 seconds, with zero failures, errors or
skips in the parsed JUnit report. This separately exercises all ten default
opt-in skips and the module CLI partition/input-byte preservation fixture.
It does not rewrite the full-suite skips as passes. JUnit SHA-256:
`aab60dc79e037704ad75caa7aca92699aef85570803b40f9f3efc86ce0d11b72`.

The native checks overlapped the unit rerun on the shared host. Durations
are regression-test observations, not controlled comparative benchmarks.
Unrelated sample-output changes were preserved. These checks do not prove
full native workflow portability, biological accuracy, corrected high-CPM
admission, resource comparability or publication readiness.
