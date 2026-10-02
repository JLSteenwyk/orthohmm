# Independent Linux Installation Verified

The new [CPU-wheel job 110652329330](https://github.com/JLSteenwyk/orthohmm/actions/runs/36947338909/job/110652329330)
at source `5b030d0fa1f9047da1d6ab89e0089233c5e01d83` completed successfully.
Build, dependency checks, both installed CLI fixtures and artifact upload passed.
This is an actual independent Linux-host installation check, not a local
checkout result or a declaration that the entire CI workflow passed.

Downloaded artifact 11202847280 once, without forwarding API authorization to
the signed artifact host. No credentials or signed redirect URLs are retained.
The [20,518-byte archive](ci_remote_cpu_wheel_artifact_20261001.zip) contains
nine report/log files, no native wheel or complete installation. SHA256 is
`0e41bf38eadf1afe7a2e16d33d22745aa00ce0063a035af289e5f640ee22271e`.
Duplicate, traversing, absolute and symlink members were rejected before extraction.
The [receipt](ci_remote_cpu_wheel_result_20261001.json) pins the actual API
observations/archive and describes its bounded verification scope.

Inspected both uploaded result and installed verification reports. Their
42 original/staged source records and driver/verifier hashes agree with Git
at `5b030d0f`. The remote package and dependency distributions are within the
fresh venv, all 36 installed package entries match the built wheel, and all
three CPU libraries load. Runtime is Python 3.12.14, NumPy 2.2.6 and Numba 0.68.0.
The declared baseline build policy is retained; no new remote ABI/ISA-wide
portability claim is inferred from the separate local ABI probe.

Both builtin/Leiden profiles produce four groups covering all 38 fixture
proteins exactly once. Both reported partition hashes agree with the local
installation fixture:
`1115fd8193636510bbc8cc8462d1b874e2a50db662fc0a59d3552d811ffa0885`.
Remote wheel SHA256 is
`ac46a766f1d07b168568c159b1159a8e12b1fd207e0653e6af86933b29efef68`
(146,122 bytes), different from the local wheel. Binary identity is not claimed.
The retained archive supplies reports/logs, not an independently re-loadable
copy of the remote wheel or a complete environment archive.

This closes independent Linux-host verification of this development-package
fixture. It does not establish frozen scientific/full phylogenetic reproduction,
dataset-wide generalization, macOS/ARM portability, accuracy, controlled
efficiency, dependency/OS closure, rights clearance or a versioned public release.
The existing macOS/full test matrix is separate and remains to be observed on
run 36947338909; prior failures remain failed, including matched-graph rendering.

The evidence-only follow-up commit uses the documented
[GitHub skip annotation](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/skip-workflow-runs)
so unchanged code is not needlessly re-executed. This does not mark its checks
passed, cancel the running source-5b030d0f workflow or change any test admission.
All implementation, workflow, tests and scientific/package files remain unchanged.
Controlled timing stays deferred, without a host poll/question or DGX access.
