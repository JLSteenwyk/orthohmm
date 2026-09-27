# Full Integrated OrthoBench Result: Job 22337

Job 22337 completed with scheduler state `COMPLETED`, exit `0:0`, 32 CPUs,
128G requested memory on `bizon`, and elapsed time `02:40:05`.
The [presubmission protocol](INTEGRATED_FULL_OB_PROTOCOL_20260927.md)
and its frozen inputs were unchanged. All eight installation, inference and
scoring stages completed in one attempt without native checkpoint reuse.

The [independent admission procedure](INTEGRATED_FULL_OB_ADMISSION_20260927.md)
then completed successfully. It checked the 215 pinned records, all stage
commands/environments, private input inventories, installed wheel inventories
and payloads, and reran all four scientific readers. It independently parsed
the full 251,378-gene partitions and recomputed all 70 reference-family scores.

## Result

| Quantity | Admitted result |
| --- | ---: |
| Genes | 251,378 |
| Root orthogroups | 59,770 |
| Weighted RefOG F1 (%) | 74.10607351873405 |
| Weighted precision (%) | 81.77045380181866 |
| Weighted recall (%) | 67.75533630827974 |
| Exact reference families | 15 / 70 |
| Changed reference-family score records | 0 |
| Genes in changed groups versus historical baseline | 0 |

All 59,770 groups are identical under label-invariant comparison. The root
partition files are also byte-identical: 6,301,075 bytes, SHA256
`5e35c1de220d0d4dca2ff64f255a55ad4d1740abbcbadbf746f91cba57f9622b`.
The independent readers report 966,439 ortholog pairs, 8,681 reconciled and
45,764 bypassed candidate families, and no checkpoint hits.

The retained [machine-readable admission result](integrated_full_ob_result_22337.json)
has SHA256 `5d3a7a1d0b84926ae7872b161520d5323514d3c22316ed067b2fe294a4c12b15`.
It includes both complete score objects, partition identities and hashes of
the detailed local provenance, package and scientific audit reports.
The local execution root is
`benchmarks/work/publication_integrated_full_ob_20260927`.

## Resources And Scope

GNU time reports 2:40:02 wall time, 268,938.42 user CPU seconds, 3,492.84
system CPU seconds, and 6,164,596 KiB maximum process RSS. These measurements
cover installation, inference and scoring on a shared host. Maximum process
RSS is not aggregate simultaneous process-tree memory. Native inference's
own timer reports 9,446.752 seconds. None is controlled comparative timing.

This closes the same-host full integrated OrthoBench reproduction check, not
independent biological validation or publication readiness. The original base
Python was used; the separately reconstructed base was tested on a fixture,
not this full run. Cross-host restoration, all-method workflow coverage,
dedicated timing, remaining uncertainty analyses, data/component rights and
public archival release remain separate requirements. Existing benchmark
scores and defaults were not replaced or selected using this result.
