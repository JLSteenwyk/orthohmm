# Installed CPU Wheel Gate

Previous turn is concrete progress: strict numerical replay environment and
copied audit evidence pushed at `3a6c2d49`. Its [CI run 36945056722](https://github.com/JLSteenwyk/orthohmm/actions/runs/36945056722)
is now terminal failure: docs succeed and all five test jobs fail. Python 3.13
reports 13,343 passed, 201 failed, 47 errors and 96 skipped. Its inspected log
shows passing overhead-audit, declared-dependency, factorial and sequence-control
modules, but matched-graph plotting still fails its strict replay guard. Do not
infer that renderer's cause from the earlier local version-only diagnostic.
Sibling failure causes were not inferred or their logs redownloaded.

The [measured receipt](ci_cpu_wheel_installation_20261001.json) pins those actual
terminal observations and the bounded new installation path below. Existing
macOS jobs remain enabled; this gate is added, not substituted for them.

## Workflow

`benchmark_tools/run_ci_cpu_wheel.py` stages only tracked package inputs. Every
file must match bytes at the recorded Git commit; dirty source, symlinked paths,
unsafe paths, duplicate names and inherited compiled binaries are rejected.
Ignored checkout libraries are never copied. The clean source is built in a
separate directory without modifying the checkout.

A fresh private venv receives the declared application/test dependencies.
The existing build backend uses `ORTHOHMM_CPU_TARGET=baseline` and build PATH
`/usr/bin:/bin`; the built wheel must be the sole new wheel. Installation uses
that artifact without resolving a different project package. `pip check` and
the unchanged installed-wheel verifier then run. Isolated `-I` probes check
package/dependency locations, every installed package byte, wheel RECORD/license
coverage and three loadable CPU libraries. Native builtin/Leiden CLI fixtures
exercise both standard and high-sensitivity profiles from the installed package.

Source/input/requirement records are checked again after execution. Logs and a
result are retained for normal errors and subprocess timeouts; unsuccessful
attempts cannot silently reuse their destination. Abrupt whole-runner termination
is not established by those unit checks. There is no job retry or cleanup of
failed outputs. Frozen scientific settings and application/build code are unchanged.

A separate Ubuntu/Python 3.12 CI job runs this workflow and uploads result/log
artifacts with `if: always()`. Native tool installation is confined to the
ephemeral CI host; no shared Threadripper package or service is changed. Its
actual remote outcome still needs observation after this prepared source is pushed.

## Actual Local Installation

This new workflow was executed once against committed package source
`3a6c2d49f201698920530fc0268ce5f6fc6e5134`, staging 42 exact files. The fresh
Python 3.12.3 venv contains the imported package and dependency distributions.
All 36 installed package entries agree with wheel bytes; three CPU libraries
load. The project license and RECORD checks pass; the selected CUDA symbol
witnesses and CUDA library entry are absent. This is not transitive rights review.

Wheel: `orthohmm-0.5.0-cp312-cp312-linux_x86_64.whl`, 146,111 bytes, SHA256
`d00277a121f2b694d13b1ad20688825df16c8fc44f6c64d286bb23f5c5e62e4f`.
No wheel was publicly released or silently substituted for a benchmark executor.
A separate isolated local ABI probe observes `hmm_have_avx2() == 0`; this does
not establish cross-architecture or whole-runtime instruction-set portability.

Both installed profiles produce four groups assigning all 38 fixture proteins
exactly once. Both partition SHA256 values are
`1115fd8193636510bbc8cc8462d1b874e2a50db662fc0a59d3552d811ffa0885`.
These are functional fixtures, not biological accuracy, full phylogenetic
inference, raw QfO/OrthoBench scoring or controlled efficiency comparisons.

**53 focused cases pass in 2.05 seconds**, zero failures/errors/skips, including
21 new orchestration cases, compilation/loading, dirty/source/input guards,
failure/timeout retention and installed-partition/wheel checks. Earlier 45- and
52-case runs are superseded focused snapshots, not additional unique tests.
Structured workflow parsing and scoped whitespace pass. This is not a new
all-suite pass, complete dependency lock, macOS fix or release readiness claim.

## Reproduction

In a Linux development verification environment with declared test dependencies:

```sh
python -m benchmark_tools.run_ci_cpu_wheel \
  --root "$PWD" --output /absolute/fresh/cpu-wheel-verification
```

The destination must be new and direct; do not remove a failed attempt to retry
under the same identity. Native GCC/OpenMP and ELF inspection tools are required.
This command installs only into its own new venv, not the controller or shared
Python environment. It is a development-package verification, not a scientific
benchmark command or permission to launch the deferred timing panel.

Other CI/native/OS/external-data failures, renderer replay, independent uncertainty,
original TreeFam files, resource execution, rights/runtime closure and final
publication release remain open. No DGX, quiet-window question/contention poll,
unrelated job/service action or new scientific benchmark occurred.
