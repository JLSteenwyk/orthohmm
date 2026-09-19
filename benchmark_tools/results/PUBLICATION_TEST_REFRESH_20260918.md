# Publication Test Refresh

Tested repository revision `e13963e` after the CLI-isolation, wheel-tagging,
native-build isolation, dependency-diagnostic and failed-output cleanup
changes. No source changes were made during the test run. Unrelated existing
sample-output changes were preserved. Tests ran on the local shared host,
not on the DGX timing node.

| Check | Result | Pytest elapsed seconds |
| --- | --- | ---: |
| Full `tests/unit` suite | 5,527 passed, 9 skipped | 122.63 |
| Native module CLI integration | 1 passed | 6.80 |
| Explicit opt-in legacy runtime checks | 9 passed | 20.41 |

The nine default skips were all controlled by
`ORTHOHMM_LEGACY_BLAST_SMOKE`. They were subsequently run explicitly against
installed native dependencies using temporary fixtures. This separate run
does not erase the skips from the original unit report.

## Commands

```sh
python -m pytest -q tests/unit --junitxml=benchmarks/work/publication_unit_e13963e_20260918.xml
python -m pytest -q tests/integration/test_module_cli.py --junitxml=benchmarks/work/publication_module_cli_e13963e_20260918.xml
env ORTHOHMM_LEGACY_BLAST_SMOKE=1 python -m pytest -q \
  tests/unit/test_admit_qfo_corrected_bpo.py::test_independent_complete_recheck \
  tests/unit/test_audit_orthomcl_database.py::test_installed_formatdb_normalization_is_detected \
  tests/unit/test_orthomcl_bpo_indexes.py::test_installed_native_bpo_parity \
  tests/unit/test_orthomcl_perl_launcher.py::test_installed_perl_probe_has_no_relative_lookup \
  tests/unit/test_prepare_orthomcl_bpo_checkpoint.py::test_installed_checkpoint \
  tests/unit/test_run_qfo_corrected_blast.py::test_installed_legacy_engine_smoke \
  tests/unit/test_stage_orthomcl_native_inputs.py::test_installed_staged_native_inference \
  --junitxml=benchmarks/work/publication_legacy_optin_e13963e_20260918.xml
```

## Evidence

JUnit reports are retained under `benchmarks/work/` without modifying raw
contents. SHA-256:

| Report | SHA-256 |
| --- | --- |
| `publication_unit_e13963e_20260918.xml` | `037ee9887f909828f0b9992f3f7367a017be321c49b9c2adf9875d0f2510b1ae` |
| `publication_module_cli_e13963e_20260918.xml` | `2a6c0e86686a59846501f29e5f610384b7d991895d3cc097f0c33e3bda288654` |
| `publication_legacy_optin_e13963e_20260918.xml` | `9de688549a88853c9ae0ac556ae43fde125e4bd66b201bc34f75c54d717ae791` |

The CLI integration checks complete, duplicate-free gene partition coverage
and input-byte preservation on a small fixture. The legacy checks exercise
BLAST success and short-query failure detection, database normalization
detection, native BPO conversion/index parity, independent checkpoint
rechecking, guarded Perl lookup and staged OrthoMCL inference. The staged
fixture expects 12 groups covering 41 proteins and 379 indexed records.

These checks do not establish biological accuracy, full-dataset completion,
cross-platform portability or comparative resource performance. Other
integration tests were not rerun here. Compiler-outcome unit tests remain
mocked; the installed CLI and legacy probes execute native tools. The
scientific executors remain frozen separately and were not updated to this
revision. No full QfO result, timing inclusion rule or publication claim was
changed by these tests.
