# Publication Test Refresh

## Current-Source Refresh At c3a89a1

Tested `c3a89a15562f15d07bfeb1c9506b703f11093be0` on 19 September 2026.
The full unit suite passed **5,842 tests**, with **9 skips**, in122.01seconds;
the native module CLI integration passed **1 test** in6.89seconds. Both raw
JUnit reports contain zero failures and errors. The only tracked change in
source/tool/test/build paths during execution was the claim-checklist prose;
no executable code changed. Unrelated sample outputs remain untouched.

```sh
/home/bizon/anaconda3/bin/python -m pytest tests/unit -q \
  --junitxml=benchmarks/work/publication_unit_c3a89a1_20260919.xml
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/bizon/anaconda3/bin/python -m pytest -q tests/integration/test_module_cli.py \
  --junitxml=benchmarks/work/publication_module_cli_c3a89a1_20260919.xml
```

Raw JUnit SHA-256 values:

- Unit: `14db6434e6bf86d60dbf0b0b6b7eefbfcea1810829c8e816f7b808a5f46bb54d`
- CLI: `ac60d1ec6ea60f9bf23569e57f7f2fd82a653e4e6bc25eb3739b1a2ccf98c3b7`

The nine skipped legacy-runtime cases were then explicitly executed and all
**9 passed in20.42seconds**. This separate run preserves the original unit
report's skips rather than rewriting its outcome. It used exactly the seven
test selectors in the earlier opt-in command below, with
`ORTHOHMM_LEGACY_BLAST_SMOKE=1`, `OPENBLAS_NUM_THREADS=1`,
`OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, the absolute interpreter above,
and output `benchmarks/work/publication_legacy_optin_c3a89a1_20260919.xml`.
SHA-256: `b208c5a1526422a32764fb302d5a167f1efb2705efd0f9f8dfe9f2a88559f9bd`.
These execute installed BLAST/database/BPO/Perl/OrthoMCL fixture checks,
not the pending corrected full-data OrthoMCL analysis.

This refresh includes the new native-pressure collector, frozen overhead
launcher and exploratory inferred-tree comparison tests. It is local shared
host validation, not a controlled timing measurement or proof of biological
accuracy, DGX outcome validity, hardware portability or publication readiness.
The dedicated DGX quiet window remains intact; scientific executors and
their frozen inputs/settings have not changed.

## Current-Source Refresh At a41256b

Tested `a41256ba6ab6942b69738b51d1b0ef4fd62d2e04` after the completed
search-decision diagnostic/join, corrected-QfO score exports, sequence-control
uncertainty and figure integration. The full unit suite passed **5,669 tests**,
with **9 skips**, in 124.50 seconds. The native module CLI integration passed
**1 test** in 7.01 seconds. Raw JUnit records report zero failures or errors.

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -m pytest -q tests/unit \
  --junitxml=benchmarks/work/publication_unit_a41256b_20260918.xml
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -m pytest -q tests/integration/test_module_cli.py \
  --junitxml=benchmarks/work/publication_module_cli_a41256b_20260918.xml
```

Raw JUnit SHA-256 values:

- Unit: `74d84ebe3d4a91a995a36bc2254fab239d2754d3ac71a04753f7e3f0621e09a5`
- CLI: `bfd233071b987d4d49de1ca58f73e5692c29bf3c64c484ee6aa1b5ff4d478ca9`

The nine skips are installed-runtime opt-in checks for legacy BLAST,
database normalization, BPO/native OrthoMCL and Perl lookup. They were not
rerun at this revision; the older explicit execution below remains separate
evidence. No tracked source, benchmark-tool, unit/integration-test, setup or
pyproject changes existed before or after these runs. Unrelated sample
outputs were left untouched. Both runs used the shared local host, not DGX.
These tests do not establish full-data completion, biological accuracy,
GPU hardware correctness, cross-platform equivalence or comparative timing.
Frozen scientific executors and their predictions remain unchanged.

## Current-Source Refresh At c7bce33

Tested `c7bce33c5fedf1af22ece559b15a544119239beb` after the search-routing
fix, retained-hit witness audit and search-coverage exports. The full unit
suite passed **5,586 tests**, with **9 skips**, in 123.10 seconds. The native
module CLI integration passed **1 test** in 6.89 seconds. No tracked changes
were present in the source, benchmark tools, unit/integration tests or
`setup.py` before or after execution. Existing sample-output changes were
left untouched. These runs used the local shared host, not the DGX.

```sh
python -m pytest -q tests/unit --junitxml=benchmarks/work/publication_unit_c7bce33_20260918.xml
python -m pytest -q tests/integration/test_module_cli.py --junitxml=benchmarks/work/publication_module_cli_c7bce33_20260918.xml
```

Raw JUnit SHA-256 values:

- Unit: `8ab7b81cfb41a8a62d1507b441b4c56569d52cedbbe5220b24ebe889eba27715`
- CLI: `786ffd8272cb2d21ee1eb18f24439cdd9aa57942509fb664287b99feaedb0e6e`

The opt-in legacy checks below were not rerun at this revision. This refresh
does not establish biological accuracy, GPU execution on physical hardware,
cross-platform portability or matched-resource performance. The frozen
scientific executors and running analyses were not modified.

## Earlier Refresh At e13963e

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
