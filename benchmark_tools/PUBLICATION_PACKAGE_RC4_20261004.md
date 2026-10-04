# OrthoHMM Study Candidate 2026.10.04-rc4

This local working candidate adds the third manuscript source and its reviewed
35-page presentation, development-family inventory, candidate reproducibility
and retained resource/provenance addenda. It preserves all 137 rc3 selected
payload identities, moving the former root guide to `evidence/`. The old
selection and reader are also preserved under `history/rc3/`, retaining the
identities of all 139 rc3 indexed payloads. The old index is included there
as an additional historical record. The frozen scientific method is unchanged.
Publication readiness, complete scientific validation, redistribution clearance,
public deposition and journal submission are not certified.

## Restore And Verify

Keep archive/index digests outside the bundle. Use the copied stdlib reader
and a new destination:

```sh
python3 -I -S -B /trusted/bundle_publication_package.py restore /candidate.tar.gz \
  --archive-sha256 ARCHIVE_SHA256 --manifest-sha256 PACKAGE_INDEX_SHA256 \
  --output /fresh/package
python3 -I -S -B /fresh/package/bundle_publication_package.py verify /fresh/package \
  --manifest-sha256 PACKAGE_INDEX_SHA256
```

No native tool starts. Verification checks payload identities and inventories,
not scientific conclusions. Generate replay output outside the immutable
package; verify the package again afterward.

| Location | Scope |
| --- | --- |
| `current-review/document.pdf` | Reviewed 35-page third-revision main text and figure appendix |
| `current-review/main.pdf`, `main.html`, `*.json` | Exact main render and source-bound presentation records |
| `evidence/PUBLICATION_MAIN_TEXT_20261004_v3.md` | Current evidence-integrated source |
| `addenda-replay/` | Six compact retained reports and standalone reporting worker |
| `evidence/` | New interpretations, execution receipts and preserved rc3 assets |
| `history/rc3/` | Original rc3 index, selection and reader identities |
| `review/`, `figure-replay/`, `statistics/`, `archives/`, `anchors/` | Unchanged inherited routes and older review |

The inherited `review/document.pdf` remains its original older snapshot,
not the current manuscript. Old guides' pending-render statements are
historical; the current review's execution receipt governs its actual status.
External file links inside PDFs/HTML and absolute paths in old records are
historical references, not portable access to those files. This bundle does
not rewrite them or claim full transitive study closure.

## Addenda Reporting Replay

Only Python's standard library is required. The worker opens its six pinned
copied reports and itself, not historical paths in their metadata. Set
`PACKAGE` to the restored directory and choose an absent external output:

```sh
python3 -I -S -B "$PACKAGE/addenda-replay/replay_publication_addenda.py" \
  --inputs "$PACKAGE/addenda-replay" --output /fresh/addenda-replay
cmp /fresh/addenda-replay/family_evidence.tsv \
  "$PACKAGE/addenda-replay/family_evidence.tsv"
```

`summary.json` recomputes exposure projections for all 88 canonical families,
the 9,156-row association TSV, fixture satellite identities/counts, native
three-repeat medians/ranges, all 24 prior register rows and eight QfO secondary
means. There are 64 metric positions plus eight secondary means, not 72
independent endpoints. It distinguishes five selected stage associations from
four observations and preserves all sixteen unavailable original full costs.

The replay does not rescan the historical Git snapshot/114 local reports,
recompute raw family-count digests, run the fixture engine, read raw accepted
merge traces, validate native output partitions or remeasure resources.
Original admissions and execution receipts retain their separate scope.
The local-only reports are pinned, not bundled. Full exposure discovery needs
those exact files plus the original Git snapshot; this report replay is not a
replacement. No independent accuracy, uncertainty or causal tuning validation
follows from arithmetic agreement.

## Preserved Routes And Remaining Work

Use `evidence/PUBLICATION_PACKAGE_RC3_20261004.md` for unchanged mechanism and
earlier comparison/factorial/native-asset routes. Their successful execution
receipts remain evidence; do not rerun completed diagnostics or timing merely
because this version exists. Shared-host timing on the local Threadripper is
authorized, but contention is unknown and potentially method-dependent.
No dedicated DGX or disruption of unrelated jobs is required.

QfO/OrthoBench remain primary and development-exposed; Three Kingdoms is
supplementary. Original TreeFam-A family reconstruction, complete causal
development history, several QfO uncertainty endpoints, prespecified biological
strata, original full configuration costs and all-method native restoration
remain incomplete. Rights, public deposition and journal formatting are also
unfinished. Preserve older archives; narrow unsupported scientific claims.
