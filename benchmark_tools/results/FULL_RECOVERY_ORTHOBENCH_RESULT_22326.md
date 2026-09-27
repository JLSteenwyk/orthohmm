# Full OrthoBench Recovery Completed

Native job 22326 completed with exit 0:0 in 2:37:28 on 32 allocated CPUs.
Independent audit job 22327 completed with exit 0:0 in 2:20 on two CPUs.
The [machine-readable result](full_recovery_orthobench_result_22326.json)
links the retained execution, installation, scientific readback and timing
receipts. The prespecified [protocol](FULL_RECOVERY_ORTHOBENCH_PROTOCOL_20260926.md)
and [readback gate](FULL_RECOVERY_READBACK_PROTOCOL_20260926.md) were unchanged.

## Result

The new hash-locked recovery installation processed all 251,378 genes from
12 raw FASTA inputs through HMM search, high-sensitivity refinement,
candidate-only canonical ordering and freshly inferred phylogeny. No search,
alignment, species-tree or gene-tree checkpoint was reused.

All 59,770 root groups match the historical partition, including exact file
bytes (SHA256 `5e35c1de220d0d4dca2ff64f255a55ad4d1740abbcbadbf746f91cba57f9622b`).
There are zero changed groups or genes. All 70 reference-family score records
match exactly, with 15 exact RefOGs:

| Metric | Recovered | Historical |
| --- | ---: | ---: |
| Weighted RefOG F1 (%) | 74.10607351873405 | 74.10607351873405 |
| Precision (%) | 81.77045380181866 | 81.77045380181866 |
| Recall (%) | 67.75533630827974 | 67.75533630827974 |

The four independent readers validated structure, sequence content,
reconciliation events/pairs and hierarchy. The result contains 54,445 candidate
families, 8,681 reconciled families, 45,764 bypassed families, 966,439 native
ortholog pairs and 200 species-tree marker families. Installation audits before
and after inference are identical. Post-completion admission rechecked the
frozen plan and its records; the linked scientific report records were also
rechecked before consolidation.

Post-completion verification: all 30 focused tests in
`test_run_full_recovery_orthobench.py` and
`test_readback_full_recovery_orthobench.py` pass. No inference source was
modified for this result consolidation.

## Resources And Scope

GNU time reports 9,438 wall seconds, 268,824.98 user seconds, 3,254.68 system
seconds and 6,163,636 KiB maximum process RSS. Scheduler duration includes
launcher overhead. Maximum process RSS is not simultaneous process-tree peak
memory. These shared-host measurements are descriptive, not controlled speed
comparisons.

This closes the same-host, full-input reproduction gap for this explicit
entrypoint and recovery environment. It does not establish cross-platform
portability, independent generalization, a new accuracy gain, or general
superiority over OrthoFinder. Canonical ordering remains an explicit
experimental policy; production defaults were not changed. Earlier discrepant
runs and historical scores remain retained, not overwritten.

The publication goal remains incomplete: controlled dedicated-host scaling,
remaining QfO uncertainty/provenance requirements and release/archive work
remain open. DGX work remains deferred. The dated manuscript review v21
predates this completed run and is not a current rendering of this result.
