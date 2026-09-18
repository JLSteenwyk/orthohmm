# Corrected Proteinortho Admission Submitted

Job `21715` is queued with `afterany:21708`, verified through `scontrol`.
It uses 2 CPUs, 32 GiB RAM and a four-hour limit on bizon, with no requeue.
The validator runs even if inference fails, but rejects any execution
without the required successful final status; failure is not converted to
a valid prediction. This job does not restart inference or score outputs.

Executor: `benchmarks/work/publication_qfo_corrected_proteinortho_admission_v1`,
detached at `8def2679d97a5507678510fa7188836006f053aa`.
The launcher checks the exact revision and tracked benchmark source state.
Batch: `qfo_corrected_proteinortho_admit_batch_20260918.sh`, SHA-256
`9cec17774bc5de50db27878c25e6794c57f27ee3134ee98eb07923bb908eb9af`.
Shell syntax validation passed.

Full unit-suite validation after submission:
`/home/bizon/anaconda3/bin/python -m pytest tests/unit -q` completed with
2,712 passed and 1 skipped in 72.00 seconds. This is software validation,
not evidence that queued inference or admission has completed.

Output, only on successful validation:
`benchmarks/work/qfo_corrected_proteinortho_admission_20260918.json`.
Log: `benchmarks/work/qfo_corrected_proteinortho_admit_21715.log`.
The output did not exist at submission. No duplicate inference was launched.

## Progress Ledger

Previous turn: progress, admission implementation pushed as `8def267`.
Current turn: frozen executor created and dependency-bound admission queued.

The corrected Proteinortho search log reported 2,298 of 3,003 species-pair
comparisons (76.52%) at inspection. This is search-stage progress, not an
end-to-end completion estimate. OrthoHMM, SonicParanoid and original QfO
factorial jobs remain active; corrected BLAST remains pending resources.

Dedicated DGX array task `21656_14` completed with scheduler exit `0:0` and
elapsed `00:44:03`. Task `21656_15` is running, with 16-26 pending. Thus
15 of 27 timing tasks are scheduler-complete; none of these observations
substitutes for final native/resource/host admission. No DGX bulk reads or
new workloads were introduced during timing.

Next: inspect `21715` after it is terminal, admit its report only if all
checks pass, then freeze and execute corrected conversion/scoring. The
full publication goal remains open.
