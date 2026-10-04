# Retained Development-Family Evidence Inventory

## Frozen Scope

Use Git snapshot `f145fff1bd1e425cdd3630ae7abdbde014fd4be0`, containing
1,785 JSON result blobs under `benchmark_tools/results` and
`benchmarks/results` (144,384,742 bytes). Also include the 114 separately
frozen local-only JSON reports under `benchmark_tools/results`
(108,265,645 bytes). Their [manifest](development_family_local_manifest_20261004.json)
has SHA-256 `905f61975705fda057a6cd4c6e42bd3d2d523d793abbaaa548f4ec32bd2eecbc`.
Original local files remain unchanged/untracked; only this new manifest and
derived inventory are delivered. Do not omit local parameter sweeps simply
because they were not committed.

Canonical family identities come only from the exact original 70-RefOG
partition manifest and the retained 18-family corrected SwissTrees count
catalog; both hashes are fixed in the reporter. Preserve the original
35/35 development/validation split as recorded chronology, not an untouched
validation certificate. Treat the primary benchmarks as development-exposed.

Explicit filename-bound `refog_records` with TP/FP/FN/genes and SwissTrees
`counts_without_prior` records define scored evidence blocks. Every block
retains its file, JSON pointer, family identities, count-vector digest and
reported partition. Plans, named/list family mentions and outcome-only
contrasts remain separate; no score is inferred from an integer family count.

Historical candidate diagnostics generate `RefOG001`-style labels by subset
ordinal (`candidate_family_diagnostics.py`), not source filename. Preserve
these containers and labels as unresolved, without assigning their counts to
canonical families. Scored synthetic runtime fixtures with `family0.txt`
labels also remain unresolved rather than becoming biological reference IDs.
No diagnostic or official score is changed or recomputed.

## Execution And Verification

Prepare/push tested source before the actual selected export. Use a standard
library-only interpreter; no native tool, scheduler job, new package or
scientific raw-data re-admission is required.

```bash
/usr/bin/python3 -I -S -B benchmark_tools/inventory_development_families.py \
  --repo . --revision f145fff1bd1e425cdd3630ae7abdbde014fd4be0 \
  --local-manifest benchmark_tools/results/development_family_local_manifest_20261004.json \
  --local-manifest-sha256 905f61975705fda057a6cd4c6e42bd3d2d523d793abbaaa548f4ec32bd2eecbc \
  --output benchmark_tools/results/development_family_inventory_20261004
```

The reporter reads immutable Git blobs via the batch object reader, not dirty
worktree versions. Local-only files are checked before reading and after the
scan, with the frozen manifest checked too. Refuse existing outputs/symlinks.
No historical absolute path inside a JSON document is opened. Preserve all
snapshot/local identities, relevant recorded source/command context, explicit
score blocks, mentions, unresolved containers, an 88-family summary and
family-by-evidence TSV. Keep duplicates; a repeated count vector is not proof
of a duplicate native execution, and a different vector is not an independent
run. Native-experiment counts and causal tuning influence remain unestablished.

Portable tests cover aliases, namespaces, numeric counts, duplicate families,
schema classification, JSON pointers, dirty/untracked snapshot isolation,
explicit local inclusion/checksum changes, unavailable independence claims
and destination refusal. Independently verify the actual output without
importing the reporter: snapshot file identities, local pins, extracted
families at every pointer, all associations/TSV fields and per-family totals.

## Development History And Remaining Scope

Before any selected export, retain these unsuccessful reporting probes:
isolated `-c` import without a repository path; the initial assumption that
all `refog_records` were official scores; and the following assumption that
all explicit count labels belonged to the biological catalog. The latter
failures exposed ordinal diagnostics and synthetic fixtures respectively.
The revised in-memory full scan passes before source preparation. No failed
scientific run or contention-based retry occurs; no preliminary result export
is selected or created.

This inventories retained family evaluation evidence, not complete causal
model-selection history or all historical commits/local native directories.
Original TreeFam-A family labels remain unavailable. Do not invent family
maps for pooled TreeFam relations, VGNC, GO/EC/FAS, BUSCO or simulations.
YGOB's admitted pre-score freeze/novel-taxon and development-homology screen
remain separately documented; no family-disjoint claim follows here. No new
method tuning, exclusion, uncertainty estimate or biological accuracy is
authorized. Existing results/defaults, native timing panel, main PDF and rc3
archive remain unchanged. Broader scientific and distribution work continues.
