# Corrected OrthoMCL BLAST Interim Observation

This is an incomplete running-search observation, not search admission,
final failure analysis or a new benchmark result. No search setting, job,
sequence or dependency was changed.

At 07:08:08 and 07:08:38 UTC on September20, controller queries confirmed
job21713 RUNNING on bizon, with zero restarts and requeue disabled. The
allocation is180CPU/900GiB and its actual controller time limit is14days.
The second observation reports elapsed18:59:15. The same native blastall
PID693830 remained present with unchanged reported start time.

The partial hit table grew from7,702,831,550 to7,707,105,684bytes:
4,274,134additional bytes over approximately30seconds. The final renamed
`all.blast` did not exist at either observation. This establishes output
activity in the observed window, not a completed-query fraction, ETA,
complete output validity or controlled resource performance. Retained `ps`
CPU percentages are process-lifetime averages, not instantaneous core use;
RSS fields are not a validated process-tree memory peak.

## Observed Diagnostics

An immutable13,762-byte log snapshot has SHA-256
`06ce2eef819902c520b50a9010f96187c90864ffe2e8ee69d018d715b58d4d7e`.
The existing strict diagnostic parser recognizes all125lines:

| Diagnostic category | Lines |
| --- | ---: |
| Setup failure | 10 |
| Karlin-Altschul statistics failure | 8 |
| Too-short query failure | 2 |
| Selenocysteine replacement | 105 |

There are87distinct diagnostic-bearing query identifiers and10distinct
failed query identifiers. Setup and underlying error lines refer to the
same failures; do not count20failed queries. Selenocysteine replacement
messages are not automatically query failures.

All10currently observed failed identifiers occur in the historical53-query
failure inventory. This does not establish equal sequence bytes, unchanged
consequences, successful processing of the other historical failures or a
lower corrected-release failure rate. These running-log counts are lower
bounds, not final counts.

The [machine-readable observation](qfo_corrected_blast_interim_20260920.json)
retains both raw controller/process responses, file size/mtime observations,
all parsed diagnostics, shared identifiers, parser/historical-audit source
hashes and the local raw-log snapshot identity. The original running log
and hit table were not modified or rescanned in full.

## Downstream Checks

Live scheduler inspection confirms the existing chain remains queued:
21713search ->21746admission ->21748BPO ->21749admission ->21750OrthoMCL
->21752admission ->21753pairs ->21754scoring ->21755admission.
Pending dependencies are unfulfilled, not reported as permanently failed.
FastOMA21740, parameter21932 and CPM21956 remain pending with their existing
downstream validation chains intact.

Reviewed `admit_qfo_corrected_blast.py` and
`audit_orthomcl_search_table.py`: final admission requires terminal successful
execution, database sequence parity and the complete renamed hit table;
logged failed queries and incoming/outgoing/self-hit coverage are retained.
Search-evidence admission explicitly does not claim success for every query.
Final grouping/reference impact and scores remain separate. The diagnostic,
table and admission tests pass:71tests in0.45seconds. No production admission
was run against partial output and no interim accuracy score was computed.

The scientific next step is to let the unchanged search finish, inspect its
complete failure/coverage audit, and retain those limitations in the final
comparison. Missing terminal evidence never authorizes a replacement run.
