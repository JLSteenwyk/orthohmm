# Corrected OrthoFinder Scoring

The existing corrected-comparator scorer and independent score validator now
support both `orthofinder_full` and `orthofinder_sequence_only`, in addition
to Proteinortho and SonicParanoid. Existing saved scores and executors are
unchanged. This is workflow preparation, not corrected OrthoFinder accuracy.

## Frozen Inputs And Semantics

Both methods require successful two-CPU conversion and exact source from clean
converter executor `aa8da7800c4801684726151ad249da7c82a3b88d`. Full/native
predictions and MCL diagnostic cliques have distinct method identities, participants,
semantics checks and result directories. Substituting root-HOG cliques or
historical participants is rejected. Reference filtering must retain every
pair without changing its bytes or hash.

The assessment environment remains SHA-256
`e86545fd04cb644ed642cf2a31b4993fb2225ca4c05743cd03e661db129e82bc`.
The same command builder evaluates GO, EC, VGNC, SwissTrees, TreeFam-A and FAS
for event year 2020. No parameter tuning, reference changes, resumed scoring,
or original-release score reuse is introduced. Short isolated work directories
are `qc_of` and `qc_om`; prior `qc_p` and `qc_s` paths remain unchanged.

The OrthoFinder scoring executor is revision
`5c34f8baad47a9895659b43696e6887aa29dbaaf`, under
`benchmarks/work/publication_qfo_corrected_of_assessment_v1`. The independent
validator retains the original `74afad5` executor binding for the two already
completed comparators, rather than retroactively relabeling their executions.

## Independent Admission

Successful scoring writes an unadmitted result. The separate validator requires
terminal eight-CPU/bizon success, reconstructs the expected preflight and full
command, verifies frozen conversion and scorer executors, rehashes their inputs
and the entire scoring output inventory, checks the native task trace, and
validates all six endpoint records. A missing/changed file or inconsistent
participant prevents admission. Admissions are new files, never overwrites.

The score-table exporter accepts each new method only with its correct
conversion status and semantics. The existing corrected SwissTrees raw-count
auditor can then use its admitted execution inventory. No table row is filled
until actual score admission succeeds; pending methods are not zero.

## Tests And Remaining Work

The 80-test focused suite covers old/new method identity and workspace
separation, wrong semantics, pending jobs, executor revision checks, scorer
failure preservation, four-method admission orchestration, output tampering,
score export and corrected raw-count validation. Orchestration tests mock
scheduler/runtime and endpoint execution, while actual record/hash handling
runs normally; they do not constitute biological endpoint execution.

Scoring tasks request eight CPUs/64 GiB/24 hours each, sequentially on bizon.
Independent admission uses two CPUs/64 GiB/four hours, also sequentially.
Corrected inference and conversion remain upstream dependencies. Paired
uncertainty still requires the complete prespecified eight-method panel;
the six-metric mean remains a project-defined secondary summary. No matched
dedicated timing claim follows from this shared-host scoring workflow.
