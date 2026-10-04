# Explicit Manuscript Export Staged

Historical staging record: the actual helper/test are now applied after all
27 terminal reviews. See [integration](PUBLICATION_EXPORT_INTEGRATION_20261004.md)
for current source validation and remaining actual-archive work.

The live review exporter binds `PUBLICATION_MAIN_TEXT_20260927.md`. Supplying
new render/print/review receipts alone cannot export the differently named
3 October or final manuscript: the historical-root and entrypoint checks
still require that old filename. The handoff builder also retains its older
nine-page review. These are packaging requirements, not reasons to alter
the frozen timing panel or overwrite historical manuscripts.

The [unapplied patch](publication_review_main_staged_20261004.patch) adds
explicit `--main-text`/`main_text` selection with review schema3. The selected
committed Markdown identity must match the rendered source and the recorded
entrypoint; selected render/print/review identities remain chained. Schemas1/2
keep their defaults and cannot reinterpret a new manuscript selection.

All70 legacy/new cases pass in6.50s, actual JUnit time6.458s. New fixtures
export an explicitly named main text, ignore dirty working-tree prose,
relocate and verify with the copied CLI under `-I -S -B` after deleting the
fixture checkout and removing Git from PATH. CLI build and invalid selections,
schema reinterpretation and missing chain evidence are also tested. These
are synthetic review fixtures, not the final study manuscript/archive or
new scientific/native validation. The [receipt](publication_review_main_staged_validation_20261004.json)
retains exact source/test/JUnit/patch identities and this limited scope.

Patch-check passes without application to the live project. Applying only
in a separate readback directory reproduces both tested files byte-for-byte.
The live review helper remains15,293 bytes,SHA-256
`f102f5f70e2c11d78e3ddfa0e0778c21372eebc4c889bb6347154ae341b31de2`.
Patch is16,624 bytes,SHA-256
`06a8654a7241f2c846e88188eb68b149edad5e44dec5cbc34309a437469b16b5`.

After all27 frozen identities are independently terminal-reviewed, recheck
the base, integrate the actual review helper/new unit file and commit them
alongside the separate source-helper export repair. Reconcile the final
manuscript, render and visually inspect its whole PDF, then export using its
exact committed receipt/revision/main-text identities. Update the final
handoff coordinator to select that review and include the committed result
helpers; its historical fixed-page candidate remains unchanged for now.
Actual final exports/relocated workflow execution and versioned archival
packaging remain required. Do not present these fixture tests as that evidence.
