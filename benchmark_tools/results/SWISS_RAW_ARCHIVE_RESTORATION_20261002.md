# Private SwissTrees Raw Archive Restoration

Extend the [raw-source relocation](SWISS_RAW_SOURCE_RELOCATION_20261002.md)
and [explicit pytest handoff](SWISS_RAW_PYTEST_HANDOFF_20261002.md) with an
input-only, bounded archive/restore helper. Both actual private panels restore
and supply the retained raw-export regressions without original-checkout reads.
This is not the full study release, public data provisioning, a rights grant,
native biological re-admission or controlled timing.

## Format And Gates

`benchmark_tools/archive_swiss_raw_sources.py` uses standard-library streaming
gzip/USTAR with deterministic metadata. It packages the unchanged pinned
`bindings.json` first and deduplicated `inputs/<sha256>` regular files. Complete
ordered original record occurrences, including duplicates, remain in the
binding; source manifests, feature/statistic identities and raw admission
checks are not rewritten or substituted.

Archiving verifies source inventory before/after and hashes the actual bytes
copied into each member. Restoration requires independently retained archive
and binding SHA256 values, validates every member's size/hash, manually writes
only allowlisted regular paths, checks the complete restored inventory and
rehashes the archive afterward. No `extractall`, symlink/hardlink extraction,
overwrite, silent retry or expected-hash recomputation. Partial failed attempts
remain for inspection. Paths must be direct and absolute; archives cannot sit
inside their input directories.

Limits: 1 GiB compressed bytes, logical payload and actual decompressed stream,
including tar headers/padding and trailing gzip streams; 2 MiB binding;
2,000 original record occurrences. Full gzip consumption checks its footer/CRC,
and actual decompressed footprint must match the canonical USTAR record size.
SHA256 pinning is integrity validation, not a digital signature. These checks
are not OS containment or a filesystem race-security certification.

## Retained Private Artifacts

| Panel | Original occurrences | Content files | Archive bytes | Payload bytes, including binding | Actual decompressed bytes |
| --- | ---: | ---: | ---: | ---: | ---: |
| Duplication | 11 | 11 | 90,734,454 | 93,703,745 | 93,716,480 |
| Fragment | 1,774 | 1,773 | 304,827,203 | 607,448,997 | 608,747,520 |

Archives remain ignored/local at
`benchmarks/work/swiss_raw_sources_{duplication,fragment}_20261002.tar.gz`.
Archive SHA256 values:

```text
duplication 3feb5037538dc8f6ea24f45806f5c0c7ea38fa9d37d48f63432a4359e1a7fba7
fragment    a0da76b9a2035197de15050017934ee7c4d7a0b406998cbb6d71914c49de03eb
```

Binding SHA256 values remain unchanged from the relocation workflow:

```text
duplication 2009be32157feddd2d0f6f4191bcbf5e0d33416b24e9a50e402f6aa023df8e58
fragment    0319644fbb580d87ab123676a7220f89f053cac89e75a5674c4295dd4d3f72f1
```

Only parties with legitimate access should transfer these files. No raw data,
native image or archive is committed/uploaded; redistribution stays uncleared.
The archives do not contain the code/runtime or all study dependencies. For
example, with the helper checkout separately available and a new destination:

```bash
python -m benchmark_tools.archive_swiss_raw_sources restore \
  --archive /received/duplication.tar.gz \
  --archive-sha256 3feb5037538dc8f6ea24f45806f5c0c7ea38fa9d37d48f63432a4359e1a7fba7 \
  --binding-sha256 2009be32157feddd2d0f6f4191bcbf5e0d33416b24e9a50e402f6aa023df8e58 \
  --output /restored/duplication-inputs
```

Then pass each restored binding and its independent digest to the explicit
pytest options documented in the previous handoff. The `archive` subcommand
accepts `--bindings`, `--binding-sha256` and a fresh external `--output`.
Do not regenerate these completed archives merely because the goal resumes.

## Actual Verification

Add 33 cases covering deterministic restoration after deleting originals,
duplicate occurrence preservation, wrong independent hashes, no overwrite,
unexpected paths/types/PAX/duplicates/sizes/bytes, gzip truncation/CRC/trailing
padding, compressed/logical/decompressed budgets and persistent/transient source
mutation. The initial 32-case panel passes. A new transient-copy control fails
before the copied-byte hash guard: source inventories match before/after, but
the archived bytes were temporarily different. Preserve that failure and the
pre-guard source snapshot. After the guard, two intermediate panels fail only
because the persistent-mutation assertion expects the old, later error message;
correct that assertion to the earlier rejection without weakening a guard.
Retain all three failed receipts. The final overlapping **233-case configured
panel passes in 12.04s**, zero errors/failures/skips.

Copy 1,827 selected source/data files (12,957,371 bytes) and both completed
archives into fresh external staging. A fresh isolated **system Python 3.12.3**
child restores both archives using only standard-library code, loads four
staged project origins and verifies all 12/1,774 members. Another fresh
**Python 3.10.13** child executes all **27 affected regression/options cases
in 6.25s**, including the four real raw-export regressions with unchanged
assertions, zero errors/failures/skips and 30 staged origins.

Each child blocks original-checkout, `/proc` and `/sys` canaries (three blocked),
has zero subsequent forbidden opens and forbids child subprocesses. Archive
and binding digests are supplied independently, not derived from received
files. Temporary staging is removed. These 27 copied cases overlap the local
233; do not claim all 233 ran in the copied tree. No native extraction,
annotation admission, inference, scoring or bootstrap is rerun. All restore
results explicitly retain `native_annotation_admission_rerun=false`,
`redistribution_authorized=false` and `publication_ready=false`.
[Machine-readable identities and receipts](swiss_raw_archive_restoration_20261002.json)
retain the exact scopes, failures and byte pins.

## Prior CI And Remaining Work

Observe source-6c55d3fd run 36997142693 at 11:07:57 UTC: docs/Linux/wheel succeed;
all five macOS test jobs fail. No restart. Actual Python 3.13 job 110806353451
log, downloaded once at 11:10:30 UTC, verifies the checkout SHA and **14,037
passes/four failures/119 skips/30 warnings in 427.61s**. All 13 new pytest option
cases pass. The four default raw exports still fail at original paths, because
public CI has neither private data nor binding options. This confirms the prior
option patch, not the new archive helper or successful full CI. Do not infer
sibling test counts from job conclusions.

Scientific engines/settings/counts/scores, two exporter helpers, raw manifests
and historical archives remain unchanged. Prospective handoffs must include
the new helper identity; old receipts are not repinned. Timing stays deferred
without further contention questions/polls, DGX or unrelated job/service
actions. Other-QfO uncertainty, controlled resource evidence, raw/transitive
rights, complete versioned release and archival deposition remain open. The
original publication goal remains active and incomplete.
