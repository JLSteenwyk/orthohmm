# Recovery Installation Advisory Range Review

A fresh authenticated, read-only [GitHub repository snapshot](recovery_dependency_alerts_20260927.json)
contains 11 open alerts, all against the retained historical
`publication_cpu_wheel_requirements_20260919.txt`, for pip and setuptools.
No alert was dismissed and no historical lock was modified.

The [new recovery range check](recovery_advisory_ranges_20260927.json) verifies
the retained recovery install-report, lock and prior byte-audit identities,
then queries the live recovery venv using its isolated interpreter. All eleven
distribution names/versions match the installation report exactly, with no
missing, additional or duplicated distribution records. The existing range
evaluator finds zero affected installed versions across the 11 reported
advisories and no alerted package absent from the environment. Recovery uses
pip 26.2.1 and setuptools 83.0.0; this result does not change the historical
installer's recorded exposure.

This checks only repository-reported ranges. It is not an all-package advisory
search, current wheel-byte re-audit, reachability/exploitability assessment,
native-library/OS review or complete build-chain clearance. Repository alerts
remain open. The full OrthoBench run and its environment were not modified.

Nine focused tests pass, covering the existing range evaluator and normalized
inventory matching plus missing, changed, duplicate and extra distributions.
An initial test command named a nonexistent test file and ran no tests; the
correct existing `test_dependency_audit.py` and new inventory tests were then
run successfully. No installation or native inference was retried.

```bash
python -m benchmark_tools.snapshot_dependency_alerts --git-credential \
  --output /absolute/new/alerts.json
python -m benchmark_tools.audit_recovery_advisories --repo . \
  --snapshot /absolute/new/alerts.json --output /absolute/new/ranges.json
```

The snapshot collector uses the existing authorized repository credentials
without including credentials in its output. These commands only read remote
alerts and local metadata; neither upgrades packages nor changes alerts.
