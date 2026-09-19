# Patched Installation Tools

The authenticated read-only [GitHub snapshot](installer_alerts_before_20260919.json)
reports11open alerts, all against pip23.0.1/setuptools65.5.0 in
`publication_cpu_wheel_requirements_20260919.txt`. These versions were
retained from the clean-venv bootstrap, not selected as a secure release
toolchain. The earlier docs dependency fixes are not the source of these
alerts. The installed `gh` is not GitHub CLI; the existing REST collector
was used without printing credentials or dismissing alerts.

The historical lock and installation evidence remain byte-for-byte
unchanged. **Do not use that lock for new installations.** For the separate
baseline CPU development artifact, use the
[patched hash lock](publication_baseline_patched_requirements_20260919.txt).
It pins pip26.2.1 and setuptools83.0.0, while preserving the eight scientific
dependency wheel versions and the previously tested baseline OrthoHMM wheel.

The [pip advisory](https://github.com/advisories/GHSA-qwm4-qh6w-59xr)
identifies26.2.0as patched; the
[setuptools advisory](https://github.com/pypa/setuptools/security/advisories/GHSA-h35f-9h28-mq5c)
identifies83.0.0. The existing range evaluator checks all11retained advisory
ranges against the actual install report: **zero affected installed
versions**. This is not a comprehensive vulnerability or exploitability
audit. Original alerts are retained, not claimed closed.

## Validation

Downloaded the two replacement wheels from official PyPI, combined them
with the retained eight dependency wheels and unchanged baseline wheel in
`benchmarks/work/publication_baseline_patched_install/wheelhouse`, and
created a fresh Python3.10.13 venv. Installed with:

```sh
"$ROOT/venv/bin/python" -I -m pip --isolated --disable-pip-version-check install \
  --no-index --find-links "$ROOT/wheelhouse" --require-hashes \
  --only-binary=:all: --force-reinstall --no-cache-dir \
  --report "$ROOT/install_report.json" --log "$ROOT/install.log" \
  -r benchmark_tools/results/publication_baseline_patched_requirements_20260919.txt
```

ROOT is the absolute path to the directory above. All11installer URLs were
local file URLs; every file size/hash was recomputed against the installer
report, and the wheelhouse contained exactly those11wheels. The initial
venv bootstrap still used pip23.0.1 to install these checked local wheels;
this is not a claim that every historical build/bootstrap tool was patched.
No global or frozen scientific environment was upgraded.

Pip check passes. The unchanged isolated installed verifier confirms
pip26.2.1/setuptools83.0.0, correct package bytes/import isolation and all
three native libraries. Standard/high-sensitivity fixtures each produce
4groups/38genes, exact coverage and unchanged inputs. Their partition
hashes equal the earlier baseline-wheel fixtures.

Retained evidence:

- [Range and local-wheel audit](patched_installer_audit_20260919.json):
  `683db1e54a55d196b71f09b2c2e3e26d244fc2b3e27f6917fa1f660f715b3b19`.
- [Installed smoke verification](publication_baseline_patched_verification_20260919.json):
  `5fd3c6cd59fa0e3a56a5a9f17d69c67f2a0e55afd6abda653ae87035baa8cb64`.
- Patched requirements:
  `77054a100352c5e676e660f0846c0ac5f1ad972786d2760f253bf87588ba1847`.

These are new development-installation results, not a retrospective change
to benchmark provenance, a complete secure build chain, a portable public
release, or closure of publication requirements. Historical installer
environments and the old lock retain their known advisory exposure.
