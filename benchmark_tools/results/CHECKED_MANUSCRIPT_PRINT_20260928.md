# Checked Manuscript Printing

The prior review exposed an unsafe orchestration pattern: printing continued
after HTML rendering refused an occupied output. Use
`benchmark_tools.print_manuscript_review` instead of an unconditional browser
command. It requires a successful renderer receipt, verifies its direct input
hashes before and after printing, and exclusively creates a fresh attempt
directory. Existing directories, files and dangling symlinks are rejected.
The browser writes only within the new attempt. A failed attempt retains
`print.json` with failure status; missing, empty or invalid PDFs are not admitted.

Example, from the repository root after a successful fresh HTML render:

```bash
python -B -m benchmark_tools.print_manuscript_review \
  --assets benchmark_tools/results/publication_main_render_20260928_v12.json \
  --output /tmp/orthohmm-manuscript-new-print-attempt \
  --browser /usr/bin/google-chrome
```

`--no-sandbox` is an explicit opt-in, used for this local trusted-manuscript
probe, not the default. Do not use it for untrusted HTML. The linked mutable
progress ledger has changed since the retained render, so a fresh render receipt
is needed to run this example successfully against current sources.

## Validation

Twenty-five printer, renderer and PDF-review tests pass. Nine printer tests
cover success, overwrite refusal, stale HTML, unsuccessful/missing receipts,
dangling output symlinks, failed browser exit, missing PDF, mid-print source
drift and timeout. A real Chrome invocation using a fresh v12 receipt produced
a six-page [PDF](publication_main_print_20260928_v12/document.pdf) with a
[print receipt](publication_main_print_20260928_v12/print.json). Repeating the
same invocation was rejected before browser launch. All six pages' extracted
text matches the visually reviewed v11 PDF exactly; this is text parity, not
a new visual review. No scientific source or benchmark settings changed.

The printer records the launcher identity, not a complete browser dependency
closure. Transitive links, network assets, visual layout and scientific claims
remain outside this check. Browser profile/cache files remain local and are
not release artifacts. Controlled timing and publication readiness remain open.
