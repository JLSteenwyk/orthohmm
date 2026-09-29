# Candidate Wheel Dependency Metadata

Offline checks bind each wheel to its locked hash, filename and top-level
METADATA name/version, then evaluate Requires-Python and Requires-Dist with
packaging's requirement, version and environment-marker APIs. The target is
the recorded executing CPython 3.10.13 Linux x86_64 environment, without extras.

| Lock | Wheels | Declared requirements | Active satisfied | Inactive markers | Failures |
|---|---:|---:|---:|---:|---:|
| Recovery | 11 | 91 | 10 | 81 | 0 |
| Reader | 5 | 52 | 2 | 50 | 0 |

[Recovery report](recovery_wheel_dependencies_20260928.json) and
[reader report](reader_wheel_dependencies_20260928.json) retain every requirement,
marker status, target version, wheel identity and Python constraint. All 20
lock/source/wheel record entries were rechecked. Thirty-one focused tests
pass, covering hash/identity mismatches, duplicate fields, missing/unsatisfied
requirements, Python constraints and inactive or unresolved extras/URLs.

Initial reads rejected pip's nested vendored METADATA files. Restricting the
reader to top-level distribution metadata fixed that error; a regression
fixture includes nested metadata. No failed attempt wrote a report or changed
an environment. A summary command had a syntax error; its corrected version
completed the identity recheck without rerunning the audits.

This verifies declared runtime dependency consistency only. It does not prove
upstream declarations complete, inspect vendored dependencies, validate wheel
tags or native loading, establish OS closure, assess security or reproduce
scientific results. The local OrthoHMM wheel's missing public PyPI endpoint
remains unresolved by the separate advisory audit. No installation, upgrade,
native inference or timing run occurred. Publication readiness is unproven.

```sh
python -B -m benchmark_tools.audit_wheel_dependencies --lock benchmark_tools/results/publication_recovery_requirements_20260926.txt --wheelhouse benchmarks/work/publication_recovery_install_20260926/wheels --output /tmp/recovery-wheel-dependencies-new.json
python -B -m benchmark_tools.audit_wheel_dependencies --lock benchmark_tools/results/publication_reader_requirements_20260927_v2.txt --wheelhouse benchmarks/work/publication_integrated_full_ob_20260927/reader_wheels --output /tmp/reader-wheel-dependencies-new.json
```

Output paths must be unused. Another interpreter/platform can select different
markers. Active URL or extra-bearing requirements are explicitly unresolved,
not silently treated as satisfied.
