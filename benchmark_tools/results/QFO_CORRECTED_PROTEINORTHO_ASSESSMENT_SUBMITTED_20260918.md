# Corrected Proteinortho Assessment Submitted

Submitted job **21718** for the six frozen endpoints GO, EC, VGNC,
SwissTrees, TreeFam-A and FAS. It uses 8 CPUs, 64 GiB RAM and a 24-hour
limit on bizon, without requeue. Initial scheduler state was pending.
Submission and successful preflight do not establish completed assessment.

Executor: `benchmarks/work/publication_qfo_corrected_assessment_v1`, revision
`74afad5376b7ee11fdabfba386851fd8d3c02857`.
Runner: `run_qfo_corrected_comparator_assessment.py`.
Batch: `qfo_corrected_assessment_batch_20260918.sh`.
Participant: `qfo_corrected_proteinortho`.

Bound pair-conversion manifest SHA-256:
`c0620eb8772a983956d2f3b5cb98636bbaca2ba0c2adbcd9ee0ff0805d02342b`.
Conversion job 21717 was independently required to be COMPLETED 0:0 on
bizon with its recorded allocation. All 4,695,385 pairs are retained.
Scoring-environment SHA-256 remains
`e86545fd04cb644ed642cf2a31b4993fb2225ca4c05743cd03e661db129e82bc`.

Fresh execution records are under
`benchmarks/results/qfo_corrected_assessment_v1/proteinortho/`, native scores
under `qfo_benchmark/scoring/corrected_proteinortho/`, and Nextflow work under
`qfo_benchmark/w/qc_p/`. Short work paths preserve the existing Darwin path
length constraint. No resume flag or historical score reuse is enabled.

## Validation and Limits

The production read-only preflight succeeded, validating pinned manifests,
runtime/reference/input files, pair identity/counts and unused output paths.
Thirty focused tests passed across the runner and conversion workflow,
including scheduler mismatch, original/failed conversion rejection,
changed filtered bytes, and preserved failed versus successful-but-unadmitted
execution states. The shell launcher passes syntax validation.

The runner saves preflight, command, pair-preparation scheduler evidence,
input/runtime hashes, logs and output inventory, and checks hashes again
after execution. A successful process remains explicitly pending independent
admission. Next validate terminal Slurm status, original preflight binding,
the complete native task trace and every endpoint/participant against frozen
reference metadata. No score is yet admitted and no paired uncertainty is
implied by process completion or native error bars.

## Progress Ledger

Previous turn: progress, corrected Proteinortho native and zero-loss pair
admission retained as `0c8adae`. Current turn: tested/preflighted the corrected
assessment runner, pinned its executor and submitted job 21718.
SonicParanoid and primary OrthoHMM inference, final original factorial
reconciliation and DGX timing remain active. Corrected OrthoFinder and
legacy BLAST remain queued. The full publication goal remains open.
