# Publication Test Refresh

## Current-Source Refresh At 53f55ff

Tested `53f55ff4666eda39a51b4370a136754e95f4bf5c` on September 19, 2026
after the parameter-neighborhood workflow and gated CPM replay changes.
The full unit suite passed **6,736 tests**, with **9 skips**, in 150.35
seconds. The module CLI integration passed **1 test** in 4.07 seconds.
Explicit opt-in legacy-runtime checks passed **9 tests** in 20.94 seconds;
these separate results do not erase the original suite's skips.

```sh
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/bizon/anaconda3/bin/python -m pytest -q tests/unit
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/bizon/anaconda3/bin/python -m pytest -q tests/integration/test_module_cli.py \
  --junitxml=benchmarks/work/publication_module_cli_53f55ff_20260919.xml
```

The opt-in command used the seven exact selectors in the earlier explicit
legacy command below, `ORTHOHMM_LEGACY_BLAST_SMOKE=1`, all three numerical
thread limits set to 1, the absolute interpreter above, and JUnit output
`benchmarks/work/publication_legacy_optin_53f55ff_20260919.xml`.

The unit result is retained here from the completed command output, not a
JUnit artifact: this unit invocation did not request XML. Both CLI/legacy
JUnit files were parsed and have zero failures/errors/skips. SHA-256:

- CLI: `24c36418909cd5e11f3c3723ed855d5b85104534f9a67340e734cb59b0549ef3`.
- Legacy: `808e01efc61c9f12954b502c97d3cf5568aeb17505fa309fed91f422e75b1f05`.

Tracked `orthohmm`, `benchmark_tools`, unit/integration tests and build
metadata were unchanged after testing. Unrelated sample-output edits were
preserved. Tests ran locally with numerical thread limits, without accessing
the DGX or changing frozen scientific executors. Native fixtures do not
establish full QfO completion, admitted CPM results, biological accuracy,
cross-platform portability, controlled comparative timing or publication
readiness. No benchmark score or scientific default changed.

## Current-Source Refresh At 58373c2

Tested `58373c23c652b2974a2558cb9142027a8a7d06b8` on September 19, 2026
after the full-node control, window-scale and dual-overhead audit changes.
The full unit suite passed **6,325 tests**, with **9 skips**, in 132.82 seconds.
Native module CLI integration passed **1 test** in 4.16 seconds. Explicitly
enabled legacy-runtime cases passed **9 tests** in 20.27 seconds; the original
unit-suite skips remain recorded separately. Parsed all three JUnit files:
zero failures and errors. Tracked source/tool/test/build files are unchanged.

```sh
pytest -q tests/unit \
  --junitxml=benchmarks/work/publication_unit_58373c2_20260919.xml
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/bizon/anaconda3/bin/python -m pytest -q tests/integration/test_module_cli.py \
  --junitxml=benchmarks/work/publication_module_cli_58373c2_20260919.xml
```

The opt-in command used the seven exact selectors in the earlier explicit
legacy command below, `ORTHOHMM_LEGACY_BLAST_SMOKE=1`, the three numerical
thread limits set to 1, the absolute interpreter above and JUnit output
`benchmarks/work/publication_legacy_optin_58373c2_20260919.xml`.

Raw JUnit SHA-256:

- Unit: `ffe4682da7b4d9cde88a9a0dfbe472f815875367968dd9d54885b47d255c81a6`.
- CLI: `083fa15efd0e8091f0e854a197850fa2e68d7d308be6fd00104574a1899cac7a`.
- Legacy opt-in: `4070c48770ff82e04228d23bc39bd466161aa8f83bfdba1e59060561f936801f`.

These local regression checks did not access the DGX or change its frozen
executor. They do not establish full corrected OrthoMCL completion, biological
accuracy, comparative resource validity or publication readiness. Unrelated
sample outputs were preserved.

## Current-Source Refresh At c0a8a37

Tested `c0a8a378d197c128f5ce1044eaa777b14de516b2` on19September2026
after the dual-bracket diagnostics/auditors and corrected factorial archive
and reproduction work. The full unit suite passed **6,137 tests**, with
**9 skips**, in135.55seconds. Native module CLI integration passed **1 test**
in4.17seconds. The nine legacy-runtime opt-in cases were then explicitly
enabled and all **9 passed in20.44seconds**. These separate results preserve
the original suite's skips rather than rewriting them as passes.

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/bizon/anaconda3/bin/python -m pytest tests/unit -q \
  --junitxml=benchmarks/work/publication_unit_c0a8a37_20260919.xml
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/bizon/anaconda3/bin/python -m pytest -q tests/integration/test_module_cli.py \
  --junitxml=benchmarks/work/publication_module_cli_c0a8a37_20260919.xml
```

The opt-in command used the seven exact test selectors in the earlier
explicit legacy command below, with `ORTHOHMM_LEGACY_BLAST_SMOKE=1`, all
three numerical thread limits set to1, the absolute interpreter above, and
JUnit output `benchmarks/work/publication_legacy_optin_c0a8a37_20260919.xml`.
All three parsed JUnit reports have zero failures and errors.

Raw SHA-256:

- Unit: `c2082be6b8c9c58762d018729175fff253007afc56aafa556e0e0ce651346c56`.
- CLI: `076985ca3fa408d19a40883519de6570e8b9a189bfdf0c36096b0f73007f094b`.
- Legacy opt-in: `78171a2ca95e04d14b4af12c9c82ece6348fba6e136df89ed600fb2bcd972d45`.

Tracked source/tool/test/build files were unchanged before and after testing;
unrelated sample outputs were not reverted. Execution was local onbizon and
did not access the DGX. The installed legacy fixtures do not establish full
corrected OrthoMCL completion. Tests do not admit live diagnostics, biological
accuracy, controlled comparative timing, cross-platform portability or
publication readiness. Frozen scientific executors were not modified.

## Current-Source Refresh At ef3ba74

Tested `ef3ba74379bd251b99578efce0cec60f5cafc51c` on19September2026
after the opt-in baseline CPU build and installed-wheel verification.
The full unit suite passed **5,960 tests**, with **9 skips**, in128.83seconds.
Native module CLI integration passed **1 test** in6.91seconds. Parsed both
raw JUnit reports and confirmed zero failures/errors. Tracked executable
sources and tests remained unchanged; unrelated sample outputs were left
untouched.

```sh
/home/bizon/anaconda3/bin/python -m pytest tests/unit -q \
  --junitxml=benchmarks/work/publication_unit_ef3ba74_20260919.xml
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/bizon/anaconda3/bin/python -m pytest -q tests/integration/test_module_cli.py \
  --junitxml=benchmarks/work/publication_module_cli_ef3ba74_20260919.xml
```

Raw JUnit SHA-256 values:

- Unit: `b2a3a83431aa71dc65ea99fe0999efc7d9ff7b3f215ec6fc60c8df4a34213179`
- CLI: `f421e4167c6c61dd5277ac10da6f32c01f5cbb574b5c7bbc88d08a375ccebe02`

The nine opt-in skips remain skips in this execution. This validates the
development checkout, not full scientific reproduction, cross-host wheel
portability or the ongoing timing/accuracy runs. No DGX access occurred.

## Current-Source Refresh At a174438

Tested `a17443819695b6467783acf847cb3221a9eb54e9` on19September2026 after
pressure-panel audit integration, analysis batch orchestration and the
corrected-strata plotter. The full suite passed **5,941 tests**, with **9
skips**, in122.66seconds. Native module CLI integration passed **1 test**
in6.87seconds. Both retained JUnit files contain zero failures or errors.
No tracked source/tool/test/build changes were present before or after
execution; unrelated sample outputs were not reverted.

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/bizon/anaconda3/bin/python -m pytest tests/unit -q \
  --junitxml=benchmarks/work/publication_unit_a174438_20260919.xml
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/bizon/anaconda3/bin/python -m pytest -q tests/integration/test_module_cli.py \
  --junitxml=benchmarks/work/publication_module_cli_a174438_20260919.xml
```

Raw JUnit SHA-256 values:

- Unit: `badb8318b55745cfaffc87bd40daa20af1a0c458d400bb6dc7e763cc347fe8ec`
- CLI: `21f4455a378d1d02fef67f12376a32d74ec52459ec41ac74d5897153afbce252`

The nine legacy-runtime opt-in tests were not explicitly enabled in this
refresh; their earlier9/9execution at c3a89a1 remains separate evidence.
Testing ran locally onbizon with no DGX access. Passing tests do not admit
the ongoing timing panel or unfinished QfO analyses, establish general
portability/biological superiority, or complete publication requirements.

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
