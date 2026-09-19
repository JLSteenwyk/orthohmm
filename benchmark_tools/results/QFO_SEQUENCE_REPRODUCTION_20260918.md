# Corrected Sequence-Control Statistical Reproduction

The committed three-arm SwissTrees bootstrap reproduced exactly in the
isolated analysis environment. A second run used a fresh directory outside
the repository, `/tmp/orthohmm-qfo-sequence.NUBqR3/reproduction`, and also
matched every numerical field. Source revision:
`5a7e7b6948a75fec794afac7927a53308fb5449f`.

The exporter reads 14 committed Git blobs: the numerical engine and its
imports, a standalone reproduction runner, audited family counts, admitted
result, frozen protocol, license and dependency requirement/lock files.
It does not read mutable checkout source bytes as scientific inputs. Worker
execution uses the exported modules under Python `-I -B`, with inherited
Python and dynamic-loader path overrides removed and numerical threads set
to one. The existing isolated environment has Python 3.10.13 and NumPy 2.2.6;
the exact installed distribution inventory is retained in each report.

The worker verifies the pinned count/result/protocol hashes and runs the
unchanged 100,000-draw, seed-20260923 bootstrap. Exact comparison covers all
11 numerical fields, including both contrasts, all six endpoints, nominal
and adjusted intervals, family differences and wins/ties/losses. Metadata
and source-admission status are kept separate: the regenerated numerical
output deliberately remains unadmitted rather than claiming a new raw-source
audit. All exported file hashes are rechecked after execution.

## Evidence

- [Initial isolated run](qfo_sequence_reproduction_20260918.json), SHA-256
  `8c6aa4d9d72ee00575789135b6f1db2b9e2cef771ef1c62e7f580d242abc0852`.
- [Outside-repository run](qfo_sequence_reproduction_relocated_20260918.json),
  SHA-256 `f67e58548ca37cfd074efce9f58541ffb30b4206bc7717c124bbb23e26cc180f`.
- [Admitted analysis and interpretation](QFO_SEQUENCE_UNCERTAINTY_RESULT_20260918.md).

Both worker processes exited zero. Eighteen focused tests pass, covering
every numerical field, changed inputs, output collisions and prevention of
source-admission promotion by a numerical-only check. Historical absolute
paths within retained JSON are provenance only; the worker does not access
them. This is same-host statistical reproduction, not native inference,
pair conversion, reference construction, official scoring or cross-platform
validation. It neither resolves third-party redistribution rights nor
establishes a complete publication release.

## Reproduce

Use a fresh output directory and report filename:

```bash
python benchmark_tools/reproduce_qfo_sequence.py \
  --repo . --revision 5a7e7b6948a75fec794afac7927a53308fb5449f \
  --python benchmarks/work/swiss_analysis_env_20260917/bin/python \
  --output /tmp/orthohmm-qfo-sequence-reproduction \
  --report /tmp/orthohmm-qfo-sequence-reproduction.json
```

The analysis virtual environment must already exist; this command records
it but does not install dependencies or authenticate all OS-level libraries.
The committed lock files are exported with the workflow for reconstruction.
