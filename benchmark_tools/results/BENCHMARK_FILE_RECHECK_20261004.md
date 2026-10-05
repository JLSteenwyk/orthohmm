# Selected Benchmark Files Rechecked

The [current all-tool register](all_benchmark_provenance_20261004_v2/register.md)
inherits input/output identities from retained admissions. The new
[streaming check](../recheck_benchmark_files.py) actually reads its 472 unique
selected file paths and compares their current sizes and SHA256 hashes with
those inherited identities. All 472 match: no missing, mismatched, nonregular,
unreadable or changed-during-read outcome. Total selected bytes are
33,912,479,521. The [actual receipt](benchmark_file_recheck_20261004.json)
retains every expected/observed identity, stat bracket and row association.

| Dataset | Method Rows | Unique Selected Input Paths | Selected Output Paths |
| --- | ---: | ---: | ---: |
| OrthoBench | 8 | 0 | 8 |
| QfO | 8 | 394 | 8 |
| Three Kingdoms | 8 | 54 | 8 |

These are file-path counts, not independent repetitions, proteome counts or
numbers of validated scientific stages. Shared files are checked once and
associated with every relevant row. In particular, the register does not
select OrthoBench input records in this field: their separate retained
provenance reports are not recursively re-audited here. OrthoMCL's QfO
input dictionary includes BPO and index artifacts, not original raw BLAST
or proof of FASTA consumption. Unknown historical inputs stay unknown.

The check does not reconstruct predictions, recompute scores, validate
conversion semantics or certify historical versions, input consumption,
transitive runtimes or complete archival reproduction. Current byte agreement
does not upgrade those claims. The existing register and its explicit gaps
remain unchanged; this is a separately bounded current-file observation.

## Validation

The [independent readback](benchmark_file_recheck_readback_20261004.json)
reconstructs all 24 method/dataset associations and all 472 unique identities
without importing the producer. It checks recorded descriptor/stat stability
and independently re-hashes all 24 selected prediction outputs, totaling
3,840,707,945 bytes. It does not re-hash the large input artifacts a second
time. All 920 frozen native helpers are independently unchanged.

The producer compares path and open-descriptor identity before and after
each streamed read, rejects nonregular paths and retains failures rather
than silently dropping them. These are per-file brackets, not continuous
or cross-file atomic integrity. Hashing uses 4-MiB chunks, not whole-file
loads. Named downstream-input dictionaries and ordinary record lists are
handled separately. Conflicting identities for a shared path are refused.
The output must be fresh and cannot alias selected evidence, including
normalized parent-directory aliases.

The first test run has 42 passes and one fixture failure: deepcopy preserves
shared aliases, so mutating one supposed conflicting pin changed all of them.
Fix the fixture to replace one pin, leaving the production conflict check
unchanged. Add unreadable, missing-result and overwrite/alias controls.
All 46 joined tests then pass, zero failures/errors/skips in 0.82s, including
16 new file-check tests. Preserve both the
[initial JUnit](benchmark_file_recheck_tests_20261004.xml) and
[passing JUnit](benchmark_file_recheck_tests_20261004_v2.xml).

## Execution Scope

Both the primary check and independent output hashing run on own CPU32,
outside 22431's native affinity0..31. No unrelated job is modified. The
primary check takes 58.498961848s; this is bookkeeping duration, not a tool
performance endpoint. It adds CPU/I/O activity while 22431 is running and
must not be represented as an uncontended interval. The existing whole-run
host observer remains active. Any timing effect is unknown; do not subtract
an assumed audit overhead or claim it is slight.

Reproduce into a fresh output path from the repository root:

```bash
taskset -c 32 benchmarks/work/release_alert_refresh_20261001/venv/bin/python -B \
  -m benchmark_tools.recheck_benchmark_files \
  --register benchmark_tools/results/all_benchmark_provenance_20261004_v2/register.json \
  --sha256 48be4b216f9f3776354246894c8c2a38f790fd4f075528e88f56b230a0a3c47c \
  --output /tmp/benchmark_current_files_fresh.json
```

Original pinned local paths must be present. No raw data are newly committed.
The live native attempt, missing uncertainty and label-independent strata,
historical provenance gaps and publication distribution requirements remain
in scope; full publication readiness is not established.
