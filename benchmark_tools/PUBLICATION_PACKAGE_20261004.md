# OrthoHMM Study Candidate 2026.10.04-rc1

This local versioned working package connects the frozen scientific method,
current manuscript/figures, shared-host resource reporting and native execution
assets. It is not a new OrthoHMM software version, a deposited public release,
a completed scientific requirement audit or submission readiness. The frozen
scientific revision is `7f3a9e40dd7e79f842cc2c11fb8b548f9a802806` (0.5.0).

The outer `PACKAGE_INDEX.json` records every delivered file, including this
guide. Retain its SHA-256 and the archive SHA-256 outside the package before
transfer. Relative package paths are usable after relocation; historical
absolute paths inside scientific receipts remain provenance, not required
reads. Nested default-candidate guides are dated historical descriptions.
Use this outer guide and each component's recorded index for current routing.

## Contents

| Location | Scope |
| --- | --- |
| `review/document.pdf` | Current 32-page review: twelve-page main, four caption-guide pages, sixteen original figures |
| `review/assembly.json` | Exact preservation, link and render records; not inference replay |
| `figure-replay/` | Pinned relative selection, all source PDFs/provenance and standalone PyMuPDF assembler |
| `archives/handoff.tar.gz` | Current twelve-page direct review, frozen scientific/setup source, workflow/tests, counts and comparison table |
| `archives/resource-reporting.tar.gz` | Final 27-attempt table/figure/prose reporting and copied reader |
| `archives/native-execution-assets.tar` | Retained wheels, MAFFT/FastTree, independent readers and controller from successful full OrthoBench reproduction |
| `anchors/` | Exact child indexes/manifests retained separately from nested archives |
| `evidence/` | Current main/claim checklist, scoped reviews, execution receipts, archive inventory and acquisition records |

The asset archive is local historical material, not blanket permission to
redistribute third-party code or data. Raw benchmark sequences/references,
native outputs for all methods, a base interpreter, Conda bootstrap and OS
libraries are not included. Acquisition recipes/manifests are in the handoff.
This is not hermetic, cross-host validated or a raw-accounting timing replay.

## Verify And Restore

The copied verifier uses only standard-library Python and does not need Git,
the repository, original workstation paths or scientific packages:

```sh
python3 -I -S -B /package/bundle_publication_package.py verify /package --manifest-sha256 PACKAGE_INDEX_SHA256
```

To restore the outer archive with a separately trusted copy of this reader:

```sh
python3 -I -S -B /trusted/bundle_publication_package.py restore /candidate.tar.gz --archive-sha256 ARCHIVE_SHA256 --manifest-sha256 PACKAGE_INDEX_SHA256 --output /fresh/package
```

Outputs must not exist. The reader checks archive/member names, types, bytes
and modes before extraction and then verifies the restored exact inventory.
It never starts scientific jobs or retries failed attempts. Verification does
not certify the scientific meaning or transitive dependencies of every file.

## Presentation And Resource Replay

Open `review/document.pdf` directly. Its original pages are preserved rather
than resized; figure guide links and bookmarks are internal. Nonfigure links
retain historical file-URI locations and are not portable raw-data access.
To regenerate the presentation with Python/PyMuPDF1.27.2.3, use a fresh output
outside the immutable package:

```sh
python3 -I -B /package/figure-replay/assemble.py --root /package/figure-replay --selection /package/figure-replay/selection.json --selection-sha256 4f21c1ea895d1d2fa26e345afbeb18aa961a3df6d37a45b8985474ce33de8f0e --output /fresh/figure-review
```

For resource reporting, restore the trusted nested archive into a fresh
directory, retain the anchor from `anchors/RESOURCE_MANIFEST.json`, and follow
`evidence/FINAL_RESOURCE_REPORTING_COMPONENT_20261004.md`. Its copied reader
replays tables, PNG pixels and generated prose from saved reviewed observations.
It does not open original raw-accounting paths or rerun native inference.
The recorded reporting stack is Python3.12.3, NumPy2.2.6, Matplotlib3.10.8.

All 27 planned identities are terminal-reviewed; 25 have measurements and 24
are eligible. Six cells have three eligible repeats; three cell summaries
remain unavailable. Contention distortion is unknown and potentially
method-dependent. Do not retry excluded attempts, select fastest repeats,
subtract background overhead or infer isolated efficiency rankings.

## Source, Counts And Native Inference

Restore the trusted handoff archive into a fresh directory and verify using
its copied `bundle_publication_handoff.py` and the separately retained
`anchors/HANDOFF_INDEX.json` digest. It includes 2,018 payload files and the
actual twelve-page main. See `evidence/FINAL_PUBLICATION_HANDOFF_20261004.md`
for exact identity and scope; its included old README describes the historical
default, not this explicit selection. Count replays and their NumPy pin are
under `arithmetic/`. Their previously successful execution remains retained;
this package does not claim to have rerun them.

The native asset archive is the exact archive used by completed job22377.
Its member inventory and historical file modes/symlinks are retained in
`evidence/integrated_execution_archive_restoration_20260929.json`. Restore it
into a fresh private directory with `tar --no-same-owner --same-permissions`;
do not extract over an existing installation. Historical FastTree mode0777
is preserved for equality, not recommended deployment permissions.

Acquire the raw OrthoBench inputs at the frozen upstream commit and prepare
the separately pinned private base/interpreter using the handoff's
`source/workflow/benchmark_tools/PUBLICATION_SOURCE_COMPONENT.md` sections
"OrthoBench Input Preparation", "Private Historical Conda Bootstrap" and
"Offline Historical Base Installation". Those sections give the fixed URLs,
versions/hashes and exact acquisition/rebinding commands; no latest-version
substitution or shared-environment/service change is authorized. Its earlier
twenty-helper omission warning applies to old exports, not this schema3
helper-inclusive handoff. Actual source import/collection evidence is separate
from syntax verification and from scientific execution.

Use the historical asset archive's `workflow.py` for its native reproduction
route, not a job-specific admission script with a substituted job ID. Supply
the freshly rebound `data.json`, its digest, compatible private base Python
and installer Python, and an unused output directory:

```sh
"$BASE_PYTHON" -I -S -B "$NATIVE/workflow.py" --assets "$NATIVE/assets" --readers "$NATIVE/readers" --reader-wheels "$NATIVE/reader_wheels" --reader-lock "$NATIVE/reader_requirements.txt" --data "$REBOUND/data.json" --data-sha256 "$REBOUND_DATA_SHA256" --base-python "$BASE_PYTHON" --installer-python "$INSTALLER_PYTHON" --output "$FRESH_NATIVE_OUTPUT" --cpu 32 --timeout 86400
```

This command installs separate scientific/reader environments, runs inference,
independent readback and full-reference scoring. It is a reproduction command,
not an instruction to rerun already completed job22377. Historical execution
was Linux x86-64 on the Threadripper with reconstructed Python3.10.13;
native tool/ABI and OS compatibility must hold. The outer archive does not
establish cross-platform compatibility or independently admit a new result.
Keep stage logs, failures and native outputs; do not overwrite old attempts.

Job22377 already reproduced all59,770 groups, all70 reference score objects
and four native TSVs exactly against its baseline, without phylogeny checkpoint
reuse. The retained result supports same-host reproducibility, not a new
benchmark score or controlled timing. Historical full-data elapsed time is
diagnostic and must not be pooled with the final shared-host resource panel.

## Scientific And External Status

QfO and OrthoBench are primary and development-exposed. Three Kingdoms is
supplementary BUSCO-restricted evidence. YGOB supplies bounded frozen
separate-clade transfer, not established family-disjoint validation. Appropriate
uncertainty remains unresolved for several QfO endpoints and the custom mean;
the original TreeFam-A7 trees/mapping were not found in inspected public sources.
Negative HMM, simulation and biological-application findings remain visible.
No custom aggregate or descriptive timing establishes general superiority.

The full goal/claim audit and final submission-format review remain outstanding.
No external archive upload, study DOI, public software release or journal
submission has been executed by creating this local candidate. The version
identifies the study package only; it does not change software0.5.0 or citation
metadata for the original preprint.
