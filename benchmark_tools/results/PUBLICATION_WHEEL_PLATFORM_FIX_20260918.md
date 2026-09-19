# Native Wheel Tag And Installed CLI Check

A clean `git archive e76a248` build exposed a packaging defect: the wheel
contained four compiled shared libraries but declared `Root-Is-Purelib: true`
and `Tag: py3-none-any`. Its SHA-256 was
`2754c8ac3a1e780a21fa5e882164cf427dadc5940101b691a715739f0955a6d9`.
The [wheel specification](https://packaging.python.org/en/latest/specifications/binary-distribution-format/)
defines the platform compatibility tag; native libraries cannot support an
unqualified any-platform claim.

`setup.py` now supplies a Distribution subclass reporting native modules,
including fallback builds. This conservatively gives platform/interpreter
tags rather than attempting ABI-independent custom tags. Compilation flags,
kernel code, inference settings and frozen executors are unchanged.
The regression test checks distribution purity, wheel platform tagging and
the console-script registration. Five focused packaging/entrypoint tests pass.

## Actual Build And Installation

Local evidence directory:
`benchmarks/work/publication_package_e76a248/`.
The original source archive was retained; its staged `setup.py` was replaced
with the tested fix for the second build. Fixed setup SHA-256:
`1d290cd8064c3d261c92cf5a604727828cf872608f2d403a5cf4a84acfd0a52f`.

Commands from that directory:

```sh
python3 -m pip wheel --no-deps --no-build-isolation ./source --wheel-dir ./wheels_fixed --log ./pip-wheel-fixed.log
python3 -m venv --system-site-packages venv
./venv/bin/python -m pip install --no-deps --force-reinstall wheels_fixed/orthohmm-0.5.0-cp310-cp310-linux_x86_64.whl --log ./pip-install-fixed.log
./venv/bin/orthohmm inputs -o outputs_fixed -c 1 --search_mode builtin --clustering leiden
```

Inputs were copies of the11top-level sample FASTAs from the clean archive,
not the working tree's generated outputs. The fixed installed CLI exited0
and produced4groups covering all38input genes exactly once. An isolated
Python import resolved inside the venv's site-packages, not the source tree.
Wheel inspection confirmed `Root-Is-Purelib: false` and
`Tag: cp310-cp310-linux_x86_64`, with the four native libraries retained.
[Machine-readable verification](publication_installed_wheel_verification_20260918.json).

Fixed wheel SHA-256:
`5348900eb29a39a5b468f6b4bcf9173c79df1504996a59df12d8a2d477d65aaf`.
Verification JSON SHA-256:
`6b3dafd0361be113f686b5551d97d79f9e6a08c320b43ef3a0ce49479e368136`.
Inference log SHA-256:
`5d34b3c607d60228bd4a8ce5caba5b7ee8e63413d3e68e2bca76448a54b80f7c`.

## Limits Before Release

This is not a hermetic dependency installation or cross-machine validation:
the venv inherits existing site packages and compilation uses host-native
flags. The platform tag does not encode CPU microarchitecture requirements
or establish CUDA/OpenMP library portability. A redistributable build policy,
library audit, other-platform validation and pinned clean-environment checks
remain required. Neither wheel was uploaded or declared a release. Retained
local build/inference times are not scientific performance measurements.
The full unit suite was not rerun after this packaging-only edit; its prior
5516passed/9skipped result remains tied to the preceding test change.
