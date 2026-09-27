# OrthoMCL OrthoBench Runs: Inputs, Conversion and Timing

## Two Distinct Retained Runs

| Evidence | April 2026 | July 2025 |
|---|---:|---:|
| Native groups | 23,803 | 23,804 |
| Assigned genes | 216,950 | 216,949 |
| Input genes | 251,378 | 251,378 |
| Changed protein sequences versus frozen inputs | 0 | 177 |
| Requested BLAST threads | 32 | 8 |
| Logged start (EDT) | 2026-04-20 16:03:24 | 2025-07-25 13:27:53 |
| Logged end (EDT) | 2026-04-23 19:12:39 | 2025-07-29 20:42:55 |
| Logged interval, seconds | 270,555 | 371,702 |

The [machine-readable audit](ob_orthomcl_provenance_20260926.json) keeps these
rows separate. Parameter and native-log timestamps agree, explicitly at
UTC-04:00. Both use logged BLAST E-value 1e-5 and MCL inflation 1.5. Exact
OrthoMCL, formatdb, blastall and MCL command tokens are retained. A binary path
containing a version is not attestation of the historical executable bytes.
Logged thread requests do not measure CPU use; CPU time and peak RSS remain
unknown. These wall intervals are not controlled timing measurements.

Both retained aggregate FASTAs have the same 251,378 gene identifiers as the
12 frozen input proteomes. Species membership from `all.gg` matches by complete
gene sets, including the shortened `Canis` label. Each indexed MCL partition
converts exactly to that run's native `all_orthomcl.out`, with no foreign,
duplicate or unaccounted indexed genes. This verifies retained conversion,
not the correctness of MCL inference or all upstream BLAST/BPO processing.

## Input-Identity Failure

April's sequence strings match exactly. July has **869 deleted `*` characters
across 177 proteins**, with no other sequence change or ID change. Two altered
human proteins, ENSP00000486295 and ENSP00000487059, occur in full RefOG023
membership before low-certainty exclusions. The audit preserves every changed
gene, original/staged lengths and sequence hashes; it does not modify the FASTA.
This observed transformation does not establish who performed it or when.

The initial strict audit raised an input-identity exception before producing
a complete report. Its [failure receipt](ob_orthomcl_provenance_failed_20260926.json)
was first committed at 15e5991. That receipt initially mislabeled the run April
because the exception lacked a run label; the explicit two-run comparison
identified July, and the receipt now records the correction. Original code and
metadata remain recoverable in that commit. The expanded descriptive reader
continues the remaining conversion/log checks, writes both runs' results, and
**still exits 1** because `input_identity_all_runs` is false. This is not a
promoted successful input-identity audit.

The earlier [partition comparison](OB_ORTHOMCL_PARTITION_COMPARISON_20260926.md)
found identical complete 70-family score records despite different partitions
outside reference families. The current table's admitted July partition must
not inherit April's exact-input assertion or its 32-thread timing. The shared
score does not prove sanitization is harmless or provide a counterfactual
effect bound. No historical scores or predictions were changed.

## Reproduction and Remaining Limits

```sh
python -m benchmark_tools.audit_ob_orthomcl_provenance --repo . \
  --base /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/ORTHOMCLV1.4 \
  --output /tmp/ob_orthomcl_provenance.json
```

The output must not exist. Expected exit is 1 with a completed descriptive
report and the failed identity flag. Thirty-eight focused tests pass, covering
timezones/weekdays, years/durations, aliases, mapping/index rejection, native
partition comparisons, wrapped MCL conversion and sequence differences.
No inference was rerun. Large raw BLAST/BPO and matrix intermediates were not
transitively audited, and the shared MCL reader is not an independent second
implementation of the clustering algorithm.
