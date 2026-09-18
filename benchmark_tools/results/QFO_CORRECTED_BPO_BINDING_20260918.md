# Corrected BLAST-to-BPO Admission Binding

`prepare_qfo_corrected_bpo.py` connects the independently admitted corrected
BLAST evidence to the audited BPO checkpoint component. It is implemented
and tested but **not submitted or run on corrected production data**.

The wrapper requires completed native BLAST accounting on bizon with
180 CPUs/900 GiB, successful independent search admission, the corrected
984,137-protein scope, and admitted hashes for the exact `work/all.blast`
and `work/all.fa` paths. Repeated identical provenance records are allowed;
conflicting records, missing records and partial-output paths are rejected.
Logged sequence-specific failures are preserved, not reclassified as
success because conversion finishes.

The independent admission job must have completed with two CPUs/64 GiB
on bizon. The admission source is bound to the clean executor
`publication_qfo_corrected_blast_admission_v1`, revision
`ae7de310cf7fbcba59107f2b054fc186f206eec4`. The wrapper verifies the supplied
admission SHA-256, its source, all checked input records and both component
audits before preparing the checkpoint. It checks the resulting source HSP
and directed-pair-block counts against the admitted table audit, and
rechecks inputs/helpers/outputs afterwards. Failures are retained in a
fresh preparation directory; no implicit resume is attempted.

Success is explicitly `corrected_bpo_checkpoint_prepared_pending_admission`,
not biological accuracy or permission to skip terminal scheduler and
independent checkpoint admission. No final inference, QfO scoring or
matched timing result is produced here.

## Verification

79 focused tests passed with `ORTHOHMM_LEGACY_BLAST_SMOKE=1` across this
wrapper, the checkpoint component and source BLAST admission. Tests cover
status/allocation/scope mismatches, record conflicts, frozen-source and
admission-hash mismatches, wrong scheduler jobs, post-admission input
mutation, preparation failures, differing source counts, and preserved
query diagnostics. The checkpoint test invokes the installed native runtime;
wrapper scheduler/admission tests use mocks, not a claimed production run.

## Execution Freeze Still Required

A system-Python import check failed because the helper import chain needs
Biopython. An isolated `python -I -B` check using the current Anaconda
interpreter imported Bio and NumPy, but also loaded unrelated editable
package hooks through system site initialization. Isolated mode does not
disable those system-site hooks. No packages or active environments were
modified to hide the problem.

Create and pin a dedicated Python environment, verify its package/source
and native-library inventory and clean startup, then freeze the executor
and batch command before submitting preparation behind BLAST admission
job21746. The current wrapper records the Python binary and version but
does not pretend those alone freeze its dependencies. Production submission
remains unexecuted. Corrected BLAST21713 is still pending resources.
