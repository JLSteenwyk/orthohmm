# Current-Source Regression Verification

At development revision `b5a1921327f1d4bb867f25f2ce6b5c38b7e42986`, the
complete unit suite passed **8,314 tests with nine skips in 213.47 seconds**.
The opt-in legacy-runtime probes plus isolated native module-CLI integration
passed **123 tests with no skips in 27.44 seconds**. Matching JUnit testcase
identities confirms all nine initially skipped cases passed in the opt-in
run. The 123-test run overlaps the unit suite; do not add the two totals.

Both commands used CPU affinity 190,191 and set `OMP_NUM_THREADS`,
`OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` and `NUMEXPR_NUM_THREADS` to 2.
They ran locally while production BLAST remained active, not on the DGX.
Test durations are regression-run observations, not comparative benchmarks.
No scientific source or test file was edited during these runs. Documentation
edits add a publication entry point and label historical launch examples.

From the repository root, with that environment and affinity:

```bash
python -m pytest tests/unit -q \
  --junitxml=benchmarks/work/publication_full_units_20260920.xml

ORTHOHMM_LEGACY_BLAST_SMOKE=1 python -m pytest -q \
  tests/unit/test_admit_qfo_corrected_bpo.py \
  tests/unit/test_audit_orthomcl_database.py \
  tests/unit/test_orthomcl_perl_launcher.py \
  tests/unit/test_run_qfo_corrected_blast.py \
  tests/unit/test_orthomcl_bpo_indexes.py \
  tests/unit/test_stage_orthomcl_native_inputs.py \
  tests/unit/test_prepare_orthomcl_bpo_checkpoint.py \
  tests/integration/test_module_cli.py \
  --junitxml=benchmarks/work/publication_native_probes_20260920.xml
```

Use new JUnit destinations when repeating. The integration test operates on
temporary input copies, verifies complete unique gene partitioning and checks
that the original and copied FASTA bytes remain unchanged. Other integration
tests were not run. This is not a new baseline freeze, complete native
comparator reproduction, production scoring admission or a publication-ready
release.

[Machine-readable verification](publication_regression_verification_20260920.json)
records source revision, interpreter/platform, JUnit identities, suite totals
and all skip reasons. Raw XML remains under `benchmarks/work/`, outside Git.
JUnit's suite duration differs slightly from pytest's terminal elapsed time;
both are retained without treating them as identical measurements.

The new [reproduction guide](../PUBLICATION_REPRODUCTION.md) has 21 local
links; all resolve to files. Its two commands are transcribed from the linked
executed reproduction records, not rerun in this verification. It distinguishes
committed-statistics reproduction from native execution and archival readiness.
