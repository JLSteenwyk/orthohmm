# Explicit Review Handoff Staged

Historical staging record: the actual helper/test are now applied after all
27 terminal reviews. See [integration](PUBLICATION_EXPORT_INTEGRATION_20261004.md)
for current source validation and remaining actual-archive work.

The [unapplied coordinator patch](publication_handoff_review_staged_20261004.patch)
connects the prepared source-helper and explicit manuscript-export repairs.
Legacy handoffs keep their dated review/default inventory. Explicit selection
creates schema3 with exact manuscript/ledger commits, main-text path, three
receipt paths and the selected page count. It requests helper-inclusive source
schema3 and verifies the child workflow revision, helper flag, manuscript
entrypoint, stage paths, page count and per-file revision mapping. It does not
accept new selection fields on historical schemas or replace the old archives.

After all27 frozen native identities are terminal-reviewed and the three
prepared patches are integrated/committed, the future build route is:

```bash
python -B benchmark_tools/bundle_publication_handoff.py build \
  --repo . --revision EXACT_WORKFLOW_COMMIT --source-profile native-build \
  --review-revision EXACT_RENDERED_REVIEW_COMMIT \
  --ledger-revision EXACT_LEDGER_COMMIT \
  --main-text benchmark_tools/results/FINAL_MAIN.md \
  --render-receipt benchmark_tools/results/FINAL_RENDER.json \
  --print-receipt benchmark_tools/results/FINAL_PRINT.json \
  --review-receipt benchmark_tools/results/FINAL_PDF_REVIEW.json \
  --output /absolute/fresh/final-handoff
```

These are placeholders, not an executed final export. The manuscript must
first be reconciled against actual final results, rendered and visually
reviewed. Receipt/ledger commits must match their historical bound bytes.
Component-local source delivery does not add raw data, runtime payloads,
scientific validation or archival deposition.

All82 focused legacy/new cases pass in5.04s, actual JUnit time4.988s. Synthetic
coordinator interfaces test binding/refusal. Two integration cases use the
real staged source/review builders with native-preparation/native-build
fixtures, retain frozen fixture scientific source, include a committed result
helper and verify the relocated combined export with copied CLI under
`-I -S -B` after removing the fixture Git checkout and Git from PATH.
These are actual fixture exports, not final study rendering/native execution
or final archive validation. The initial80-case run retains one mock-fixture
failure and79 passes; repair the missing committed-file reader and add the
integration cases before the successful82-case run. Both JUnits are pinned
in the [validation receipt](publication_handoff_review_staged_validation_20261004.json).

Patch-check passes without changing live project source. Isolated application
reproduces the tested helper/unit file exactly. Live coordinator remains
17,710 bytes,SHA-256
`9be5629233687af860856d297a698b7df0705d44e208c911f1937d0687f7ef6e`.
Patch is25,375 bytes,SHA-256
`f5fe6a7718ee2632b08d22f20d9351b7125ccd58f9fdbdb319ec86ef8793a72d`.
Keep the timing recipe intact while22424 runs. Recheck/apply all three actual
helpers only after terminal review, then execute and validate the committed
final exports. Scientific claims, final manuscript/reporting/release packaging
and the full publication goal remain incomplete.
