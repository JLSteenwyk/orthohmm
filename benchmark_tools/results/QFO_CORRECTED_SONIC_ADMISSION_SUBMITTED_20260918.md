# Corrected SonicParanoid Admission Submitted

Job `21716` is scheduler-confirmed pending on `afterany:21710`. It uses
2 CPUs, 32 GiB RAM and a four-hour limit on bizon, without requeue. It
does not repeat inference or score predictions. A failed parent causes
admission failure, not acceptance of partial tables.

Frozen admission executor:
`benchmarks/work/publication_qfo_corrected_sonic_admission_v1` at
`8db0c494cc97a6097ba3f178dd79496b1d29332e`.
Batch file: `qfo_corrected_sonic_admit_batch_20260918.sh`.
Expected output: `benchmarks/work/qfo_corrected_sonic_admission_20260918.json`.
Log: `benchmarks/work/qfo_corrected_sonic_admit_21716.log`.

## Gates and Validation

The validator requires the expected completed Slurm identity/allocation,
successful native execution, unchanged frozen inference executor, exact
plan/source/command/environment bindings, identical copied inputs, and
matching pre/post/current runtime inventories. It verifies the full output
file inventory and hashes, then requires all native species-pair tables,
correct input-species ownership, valid counts and finite confidence values.
The existing relation deduplication behavior is preserved. Input, output
and source hashes are checked again after traversal.

The real command was tested against running job 21710 and rejected it at
the scheduler-completion gate. No admission report was created. Forty-four
focused tests passed across the admission, table validation and runner
modules; shell syntax validation passed. Historical tables had separately
passed the full ownership audit, but do not substitute for corrected outputs.

These gates establish representation/provenance, not biological accuracy,
independent proof of algorithm completeness or matched timing. Pair
conversion, reference filtering and assessment remain separate stages.

## Progress Ledger

Previous turn: progress, native table validator and historical audit pushed
as `de12ad2`. Current turn: tested terminal admission and queued its pinned
executor. No existing inference or unrelated job was stopped or restarted.

At final scheduler inspection, 21710, 21708, 21706_0, 21711 and 21671_3
remained running; admissions 21716, 21715 and 21712 waited for dependencies.
BLAST 21713 waited for resources, OrthoFinder 21706_1 for its array slot.
DGX task 21656_15 remained active. Next: inspect terminal admission results,
then prepare corrected conversion/scoring. The full publication goal remains
open; no corrected SonicParanoid score has been reported.
