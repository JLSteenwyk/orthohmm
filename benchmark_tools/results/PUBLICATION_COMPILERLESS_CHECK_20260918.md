# Installed No-Native-Library Check

Tested the locally installed wheel from the preceding build-isolation audit,
temporarily moving its four `.so` files into a holding directory. This was
confined to `benchmarks/work/publication_package_e76a248/venv`; no active
benchmark environment or checkout library was changed. All four libraries
were restored after the two no-native runs completed.

| Installed command mode | Exit | Observation |
| --- | --- | --- |
| Standard, native libraries absent | 0 | 38genes,4groups; orthogroup file matches prior native standard output byte-for-byte |
| High sensitivity, native libraries absent | 1 | Profile construction raises missing `pair_align.so` error |
| High sensitivity, native libraries restored | 0 | 38genes,4groups,195network edges |

Each used the same copied11FASTA fixture, `-c 1 --search_mode builtin
--clustering leiden`, separate output directories, and the installed console
script outside the checkout. High sensitivity added
`--accuracy_profile high_sensitivity`. These observations do not establish
cross-platform equivalence or biological accuracy; elapsed times are not
controlled performance evidence.

Raw logs under `benchmarks/work/publication_package_e76a248/`, SHA-256:

| Log | SHA-256 |
| --- | --- |
| installed-inference-no-native.log | 8556efb9658352d5723d32f3d3820ffebd960b8e0822e53114e3336da66e0238 |
| installed-inference-no-native-high.log | 10349ee23c6e4269aad00c7b5fcf0b92fb00e8b5af94cd9420fd834d6c976202 |
| installed-inference-native-high.log | d4588e7f38866070eeee64b76decfc0b11496c417c0542f075b403db02684483 |

## Change And Scope

Qualified setup warnings and README guidance: the Numba standard-search
fallback does not implement profile alignment. Native alignment load failure
now raises an actionable RuntimeError naming the required kernel, compiler
and OpenMP support while preserving the original OSError as its cause.
No silent fallback, disabled refinement, score change or alternative aligner
was introduced. The installed experiment precedes the message-only change;
that new diagnostic is separately covered by unit tests.

Verification:9packaging/entrypoint/native-dependency tests and14profile-
expansion/CLI-integration tests pass (23total). Cached-loader behavior is
unchanged. The full suite was not rerun in this check. Compilerless support
for the complete high-sensitivity method remains absent, now explicitly
documented rather than implied by the standard-search fallback.
