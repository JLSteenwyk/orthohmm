# Documentation AnyIO Security Update

GitHub alerts44/45 identify AnyIO in `docs/uv.lock`, not OrthoHMM's
declared inference dependencies. The retained
[pre-push snapshot](docs_anyio_alerts_before_20260918.json) records:

- [GHSA-82r6-8w77-94w6](https://github.com/agronholm/anyio/security/advisories/GHSA-82r6-8w77-94w6):
  critical TLS hostname-validation issue, patched in4.14.2.
- [GHSA-5p39-cfhj-2xmp](https://github.com/agronholm/anyio/security/advisories/GHSA-5p39-cfhj-2xmp):
  medium process-worker stderr blocking issue, patched in4.14.2.

Added `anyio>=4.14.2` to the existing docs-only uv constraints and regenerated
the lock with uv0.12.15 using `uv lock --upgrade-package anyio`. The resolver
selected4.15.1 instead of4.7.0, removed the no-longer-required sniffio entry,
and updated typing-extensions to4.16.0. Other locked versions are unchanged.
No global or frozen benchmark environment was upgraded.

## Verification

- `uv sync --locked --all-extras --dev` installed the docs environment under
  Python3.12.3.
- `uv run --locked python -m sphinx.cmd.build -W --keep-going -b html . ../benchmarks/work/docs_anyio_security_20260918/html`
  completed successfully for all7documents.
- `uv run --locked sphinx-autobuild --help` completed successfully; no preview
  server was exposed.
- Existing dependency-audit tests:4passed.
- [Lock audit](docs_anyio_lock_audit_20260918.json) checks every locked branch
  against both retained ranges:0affected alerts. This is not a universal
  vulnerability scan or proof of application exploitability.

GitHub alert closure requires a separate post-push API observation; it is
not inferred from the local lock check. Historical benchmark manifests and
running jobs remain unchanged.
