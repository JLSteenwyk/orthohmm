# Native Duplication Projection Replayed From Private Restored Inputs

[Protocol](NATIVE_DUPLICATION_PORTABLE_PROTOCOL_20261006.md) pushed atd8c70fd6;
worker/tests43bfceab before actual packaging/replay. Existing original projection,
reader and frozen scientific files unchanged. New
[replay receipt](native_duplication_portable_replay_20261006_v1.json) confirms
18 families/563 genes,21,530 raw rows,36 family records, eight projections/four
differences and original TSV values. Statistics agree at original1e-12 tolerance;
no primary JSON regeneration/primary-exporter rerun or new scientific admission.

## Private Package

Existing content-addressed stager preserves34 direct source occurrences, including
one duplicate, in33 distinct input members plus `bindings.json`. Private bounded
tar/gzip archive/restore succeeds:43,904,145 payload bytes,43,939,840 decompressed
tar bytes,34 archive members. Binding bytes/hash identical after restore.

| Local Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| `benchmarks/work/native_duplication_handoff_20261006/native_duplication_inputs.tar.gz` | 39171350 | `71b435ddfc8a23c142be64c9742233c38d0a292fd656c4149d871dc0b259bd67` |
| `staged/bindings.json` and restored binding | 13713 | `a59f9208507b38de5416aca2f8504092ebadbd67c6949fbf4db89b529aa14088` |
| `source_43bfceab.tar` in the same work directory | 194560 | `f544096d1b1f513ed47c5a884ded162f840d6d3521cb3d6b55d3ce607b2569d6` |
| Replay worker | 9152 | `24b2c16994a3333e82b5d184e3f8c4a4d3727ffa19ebaac8353205d197dad095` |
| Replay receipt | 29742 | `b99d80205d39918142d448b1b2113d9510fc982c153e998188c96a765af9c770` |
| Source identity receipt | 13680 | `ad75cfe1eed7dcffc06d1026c6901d14ec623c4b8b2b36a6700124301802edd2` |

These local archives are retained, not committed or publicly deposited. No
redistribution authorization, data-rights clearance or final release assertion.
Canonical report paths remain provenance; copied inputs are content-addressed,
not rewritten into a new scientific reference. No original-path fallback.

Git archive at `43bfceab7a54359686819b87281139ec66c43cdc` includes19 files:
worker, original primary source, four replay/relocation helpers, archive helper,
canonical projection/TSV/Markdown, three protocols/amendments, five focused test
modules and LICENSE. [Identity receipt](native_duplication_portable_source_identity_20261006.json)
checks all19 Git/tar/copied bytes and compiles12 Python files. This complete
export comparison is post-replay, not a pre/post bracket for unused fixtures/docs.
Worker separately checks its operational sources/inputs/outputs before and after
computation. No continuous/OS-level integrity or security certification.

## Real Replay And Tests

The private raw archive restores under
`/tmp/native-duplication-portable-20261006-43bfceab/restored`; source tar extracts
under its separate `source` directory. Fresh isolated scientificPython3.10
process runs exported worker with-I-S-B from the separate `/tmp` working
directory, stdlib only. Exact accession/reference-entry sets and alias handling,
original count-audit/cell bindings, failed R1 timing exclusions, identical raw
truth, family statistics, integer-bin decisions, differences, null/NA empty
bin and every serialized TSV cell check. Primary exporter, tree extraction,
native inference/scoring and bootstrap are not rerun.

The original files remain physically available on this host. No fallback is
established by source behavior and synthetic replay with originals removed,
not an OS-denial/chroot audit or cross-platform validation. The interpreter's
stdlib/runtime still comes from the retained local installation.

[Final joined JUnit](native_duplication_portable_20261006/native_duplication_portable_joined_tests_20261006_v2.xml):
176 tests/zero failures/errors/skips,2.322s XML/displayed2.36s;27 new cases.
Full synthetic replay removes originals; payload/source/binding/TSV mutation,
unknown logical refs, scope flags, count/native-cell identity, failed-timing
resources/eligibility and output refusal tested. Existing bounded archive/
relocation/independent statistic tests joined. Initial19-case receipt retained
as development history, not final tested source. No scientific runtime install.

Sanitize PYTHONPATH/PYTHONHOME/PYTHONUSERBASE/LD_PRELOAD/LD_LIBRARY_PATH/LD_AUDIT;
set PYTHONNOUSERSITE=1,PYTHONDONTWRITEBYTECODE=1,PYTHONHASHSEED=0,
OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1. Original commands in copied receipts.

| Shared-Host Postprocessing Phase | Wall Seconds | Maximum Process RSS (KiB) |
| --- | ---: | ---: |
| Stage | 0.29 | 15360 |
| Archive | 1.80 | 15360 |
| Restore | 0.59 | 18432 |
| Exported replay | 6.45 | 1655260 |

All four time receipts exit0/zero swaps. These are observed shared-host
packaging/postprocessing costs, not inference timing or isolated speed.
Competing CPU/memory-bandwidth/I/O effects are unknown and potentially
tool-dependent. Pre-launch RAM658,819,556KiB available/free swap24,988KiB;
near-full swap noted, no unrelated workloads/services altered.

To replay from a restored package, use fresh output and pinned binding:

```bash
python -I -S -B SOURCE/benchmark_tools/reproduce_native_duplication_projection.py \
  --projection SOURCE/benchmark_tools/results/native_qfo_swiss_duplication_strata_20261006_v1/report.json \
  --bindings RESTORED/bindings.json \
  --bindings-sha256 a59f9208507b38de5416aca2f8504092ebadbd67c6949fbf4db89b529aa14088 \
  --scores SOURCE/benchmark_tools/results/native_qfo_swiss_duplication_strata_20261006_v1/scores.tsv \
  --table SOURCE/benchmark_tools/results/native_qfo_swiss_duplication_strata_20261006_v1/TABLE.md \
  --output /absolute/fresh/replay.json
```

Use recorded sanitized environment. Restore with existing `archive_swiss_raw_sources`
CLI using the archive/binding hashes above. Raw inputs still require private
handoff; this is not a self-contained public data release.

## Remaining Scope

Advances executable reproduction7.4 for this component, not the complete native
factorial/full pipeline or final publication archive. Broader uncertainty,
generalization/strata/provenance/TreeFam and release limits remain. Original
22444 now COMPLETED0:0,11:42:59;22444.0 COMPLETED0:0,11:18:09; original reviewer
22445 RUNNING1:36,22450/22451/22452 pending. Wait for its actual terminal/review
gates before nextnativeidentity. No unfinished reviewer output read, original
job restarted or nextidentity released. Full publication goal active/unproven.
