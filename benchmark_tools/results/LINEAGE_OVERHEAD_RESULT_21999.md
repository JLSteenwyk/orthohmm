# Complete Lineage Collector Overhead Result

All 18 array-21999 tasks and recorder 22000 completed with exit code 0:0.
Native outputs were collected only after every task and the recorder were
terminal. Independent controller replay checked 2,413 polls with zero
observation errors and retained first-terminal records for all 18 tasks.
The complete archive contains the deployed recipe, original four-proteome
inputs, native outputs, raw observations, terminal scheduler records and logs.

The pinned `lineage_21999` audit validated all 18 tasks, with no missing or
invalid evidence, temporal-order issue or changed boot domain. All nine pairs
passed output-equivalence and minimum-duration checks. Every pair met the
prespecified signed overhead limit of 10%, and each method's three-pair
median met the 5% limit. This is a collector diagnostic on the smallest
scaling input (73,266 proteins, 20 CPUs, 96 GiB), not a scientific speed ranking.

| Native method | Pair 1 (%) | Pair 2 (%) | Pair 3 (%) | Median (%) |
| --- | ---: | ---: | ---: | ---: |
| OrthoHMM high sensitivity | +0.085529 | +0.176313 | -0.372571 | +0.085529 |
| OrthoHMM satellite_v2 | +0.675138 | -0.668083 | -0.505169 | -0.505169 |
| OrthoFinder full | +1.614059 | -0.352876 | +2.668681 | +1.614059 |

Values are `100 * (periodic_wall / boundary_wall - 1)`, not overhead-subtracted
timings. Negative values reflect run variability, not evidence that monitoring
accelerates inference. No pair was omitted or selectively rerun.

## Screening Limits

All whole-command screens passed, but periodic tasks retained 143 original
and 23 narrow interval flags. Original/narrow counts by task index were:
1: 1/0; 3: 44/6; 5: 2/1; 6: 49/10; 8: 1/0; 10: 2/1; 13: 3/1;
15: 1/0; 17: 40/4. Boundary arms explicitly retain null interval flags and
unavailable interval coverage. They cannot exclude transient competing work.
No cause is inferred from these flags or non-atomic CPU-counter differences.

The numerical overhead budget passes, but environmental validity, controlled
scientific timing admission and publication readiness remain false. Remaining
CPU discrepancies, during-read lifecycle behavior and a prospective scientific
inclusion policy still require resolution. Existing historical runs are not
retroactively promoted. No accuracy score or scientific default changed.

## Reproduction

- Full local audit: `benchmarks/work/lineage_overhead_audit_21999.json`.
- Retained audit: `dgx_lineage_overhead_audit_21999_20260919.json.gz`, SHA-256
  `7e4cf99291450c883bf69abe0d984d07737953c2792dcb191bdb5a3cbc73251e`.
- Controller capture, replay and accounting: `lineage_overhead_scheduler_21999.tar.gz`,
  SHA-256 `357db8ce0db68f0788b6ccd2c7c69031361e4a84cdf4a3504053b29e2c792b6f`.
- Compact machine-readable result: `lineage_overhead_summary_21999_20260919.json`.
- Archive: `benchmarks/work/lineage_overhead_archive_21999/`; 649,502 transferred
  regular files totaling 3,769,949,218 bytes before adding scheduler records/logs.

The [audit instructions](LINEAGE_OVERHEAD_AUDIT_READY_20260919.md) give the
complete raw-replay command and provenance gates. Regenerate the compact result:

```sh
jq --arg sha256 7e4cf99291450c883bf69abe0d984d07737953c2792dcb191bdb5a3cbc73251e \
  -f benchmark_tools/summarize_lineage_overhead.jq \
  benchmarks/work/lineage_overhead_audit_21999.json
```

Validation: 297 focused summary/provenance/orchestration/arithmetic/lineage
replay tests pass. Separate Python arithmetic recomputed all nine signed
ratios and three medians from audited arm times (agreement within 1e-12),
checked the retained audit hash, flag totals and non-admission fields.
Both compressed artifacts pass gzip integrity checks. Scoped whitespace
checks pass; pre-existing whitespace in unrelated sample outputs is untouched.
