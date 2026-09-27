# OrthoBench Scoring Archive Restored

The [build receipt](ob_scoring_archive_20260927.json) records a deterministic
local archive of the validated full-data scoring component. Two independently
written archives are byte-identical: 1,383,071 bytes, SHA256
`64f193cbedfe2f03251ae539ebf355a43bd4d6961d5297bb02cd162bf04c457b`.

Local archive:
`benchmarks/work/publication_ob_scoring_archive_20260927/scoring.tar.gz`.
Its 40 regular-file members include the frozen 31-module reader export,
scoring worker, strict input manifest, retained root-group predictions,
expected score, patched reader lock, project license and standalone restore
script. It contains no raw upstream FASTA/reference files, installed venv,
wheels, native inference tools or operating-system libraries. It is not a
public release or a redistribution-cleared archive.

## Executed Restoration

The archive and its restore script were copied/extracted into
`/tmp/orthohmm-ob-score-restore-stage-20260927`. The script was invoked there
with isolated Python, a clean environment, and no checkout imports. It
verified the externally supplied archive digest before processing members,
rejected unsafe/nonregular/duplicate members, verified every member identity,
and restored 93 separately supplied raw inputs against archived hashes.
Those inputs were the previously verified acquisition copies, not a fresh
network download. The patched reader interpreter was supplied separately.

Restoration to `/tmp/orthohmm-ob-score-restored-20260927` completed once with
no retry. The worker checks the complete staged input inventory and recomputes
the score from memberships. The [restoration receipt](ob_scoring_archive_restoration_20260927.json)
records exact equality of the full score object: 251,378 genes, 59,770 groups,
70 reference families and F1 74.10607351873405%. The entire 15,378-byte score
file is identical to the earlier portable result, SHA256
`0525f414fa03578da72fb01d6116fd9696a8a230b25b8421335c672f5c4be3a3`.

The restore process and worker descendants were file-traced; neither the
original repository prefix nor the main installation's site-packages prefix
appears. This literal-prefix check is not a filesystem sandbox. All archived
payloads were rechecked after execution. Forty focused tests pass across the
archiver/restorer, portable worker, reader exporter and acquisition verifier.

## Reproduction

Use a fresh destination and the separately installed patched reader environment.
Verify the archive digest against the pinned receipt **before extracting or
executing its restore script**. Then run the trusted extracted script:

```bash
/patched/reader/venv/bin/python -I -B /trusted/restore.py restore \
  --archive /local/scoring.tar.gz \
  --sha256 64f193cbedfe2f03251ae539ebf355a43bd4d6961d5297bb02cd162bf04c457b \
  --acquisition /separately/acquired/Open_Orthobench \
  --python /patched/reader/venv/bin/python \
  --output /fresh/restored-scoring
```

`--acquisition` contains `BENCHMARKS/Input` and `BENCHMARKS/RefOGs`.
Only pinned files are copied; input bytes must match exactly. Restore refuses
an existing output directory. The source controller can rebuild the archive
with `python -m benchmark_tools.archive_portable_ob_score build`, using the
validated score root and receipt; restoration does not require that controller,
Git history, the original admission path or the earlier temporary score root.

This establishes same-host restoration of the full scoring component, not
automatic reader installation, native inference restoration, cross-host
portability, new accuracy evidence, QfO reproduction or publication readiness.
The complete acquisition/install/inference/readback/scoring workflow, public
deposition, rights review and existing scientific/timing gates remain open.
