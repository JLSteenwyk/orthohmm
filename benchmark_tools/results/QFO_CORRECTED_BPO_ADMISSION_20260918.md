# Independent Corrected BPO Admission

Submitted **21749**, dependency `afterany:21748`, with two CPUs, 64 GiB RAM,
24-hour limit, node request bizon, no requeue. Authoritative scheduler state
is PENDING/Dependency. `afterany` allows the validator to report failure if
preparation terminates unsuccessfully; it does not admit failed preparation.

Frozen validator: `publication_qfo_corrected_bpo_admission_v1`, revision
`914c3fe6f8618c3d90452e7a5b6d7f03f939ccd2`.
Batch: `benchmark_tools/results/qfo_corrected_bpo_admit_batch_20260918.sh`.
Preparation job argument:21748; preparation executor remains
`826e963cf06ed414609f0609b96c6bfd283fc2d8`.

## Admission Checks

- Require unique terminal COMPLETED/0:0 preparation accounting, the intended
  allocation/node, matching job ID, finite ordered timestamps and pinned
  before/after Python-runtime evidence.
- Recheck frozen preparation revision/source, the exact checkpoint output
  inventory and guarded native commands, input/helper/output hashes and
  empty native stderr. No changed paths or implicit cache substitution.
- Revalidate the original BLAST admission, including its executor, scheduler,
  input identities, corrected 984,137-protein scope and query diagnostics.
- Independently rerun the full HSP-to-BPO content check and native index
  validator in a fresh directory. Recheck all byte offsets, EOF sentinel and
  inclusive query ranges, not just the summary counts.
- Verify Python and native Perl/system runtime inventories before and after;
  rehash checked artifacts and retain failures. Successful output identifies
  the admitted BPO/index files but leaves accuracy/publication flags false.

Expected result directory:
`benchmarks/work/qfo_corrected_bpo_admission_20260918/`.
Expected success status: `corrected_orthomcl_bpo_checkpoint_admitted`.
No full-data checkpoint has yet been admitted.

## Verification

86 focused tests passed with installed legacy-runtime checks enabled.
Tests cover provenance/status/allocation/timestamp/runtime/command/inventory
rejection, damaged native indexes, complete native fixture rechecks, and
admission orchestration success and retained failure paths. Scheduler
orchestration tests are mocked; actual fixture checks use installed Perl.

The dedicated interpreter independently rechecked the retained frozen native
fixture: six BPO records, three queries and seven offsets, with Python
runtime verification before and after. Retained component results:

- `orthomcl_independent_bpo_content_20260918.json`, SHA-256
  `1f931942f200ae69097521dd92b34bdd41082ada94cd61fb1542d59431a73450`.
- `orthomcl_independent_bpo_indexes_20260918.json`, SHA-256
  `1a7b4921402625e586eee9a2ab9e10e482bb5fa8efb4f32f6c24425edd0b4402`.

The frozen validator was also invoked against still-pending21748: it rejected
nonterminal scheduler evidence and created no output directory. Batch syntax
and clean executor checks passed. Raw independent fixture artifacts remain
under `benchmarks/work/orthomcl_independent_bpo_recheck_20260918/`.

These checks establish conversion/index consistency, not biological truth.
They share reviewed legacy arithmetic with the converter; the native BioPerl
parity fixture is complementary, not a full independent native production
parse. Guarded native inference, final-group validation, pair conversion,
QfO scoring and failure-impact review remain required.
