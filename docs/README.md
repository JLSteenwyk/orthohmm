# Documentation Environment

The docs build uses Python 3.12 and uv 0.12.15. The docs project requires
Python 3.10 or newer so patched HTTP and live-preview dependencies can be
resolved. This does not change OrthoHMM's package runtime requirements.

From this directory:

```sh
uv python install
uv sync --locked --all-extras --dev
uv run --locked python -m sphinx.cmd.build -W --keep-going -b html . ./_build/html
```

The constraints in `pyproject.toml` prevent the vulnerable versions identified
in the September 2026 dependency audit from being selected again. They are
minimums, not a guarantee against future advisories. Update and audit the
lockfile before deployment. Keep live preview bound to localhost and do not
build untrusted documentation or templates with credentials available.

Frozen publication benchmark environments are separate historical artifacts;
do not upgrade them in place when updating this docs lockfile.
