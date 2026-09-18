# Corrected FastOMA Scoring

## Frozen Executors And Jobs

| Stage | Revision | Job | Dependency | Resources |
| --- | --- | --- | --- | --- |
| Six-endpoint assessment | `9258bcfd3f90d63ec7f2cfb02122cd20bd7e1214` | 21744 | afterok:21742 | 8 CPUs, 64 GiB, 24 hours |
| Independent score validation | `a076a36e09f5ea98b6433d53641301011be85843` | 21745 | afterany:21744 | 2 CPUs, 64 GiB, four hours |

Both jobs are pinned to bizon with automatic requeue disabled. Detached
executors are `benchmarks/work/publication_qfo_corrected_fastoma_assessment_v1`
and `benchmarks/work/publication_qfo_corrected_fastoma_score_admission_v1`.
The preceding pair conversion is frozen at
`6616e3a7ec4c46962ded0d34ce9b4150073a2210` (job 21742).

## Admission Rules

The scorer requires completed two-CPU pair conversion, the exact frozen
conversion source, corrected participant identity, positive distinct-pair
counts, raw-row/duplicate accounting, zero reference-mapping loss and byte
identity between submitted and reference-filtered pairs. The conversion's
mapping must match the frozen assessment environment. All upstream checked
records are rehashed before and after native scoring.

FastOMA retains native phylogenetically inferred pairs and the explicitly
supplied corrected OrthoFinder species tree. No HOG-clique substitution or
claim of independent species-tree inference is made. Existing Proteinortho,
SonicParanoid and OrthoFinder scoring executors remain unchanged.

The six unchanged endpoints are GO, EC, VGNC, SwissTrees, TreeFam-A and FAS.
Workspace `qfo_benchmark/w/qc_f` is distinct from other methods. Native score
files go to `qfo_benchmark/scoring/corrected_fastoma`; execution records go to
`benchmarks/results/qfo_corrected_assessment_v1/fastoma`.

The independent validator rejects unsuccessful scoring even though scheduled
after any terminal scoring outcome. It checks scheduler identity, execution
and conversion provenance, source revision, reference mapping, exact commands,
output inventory, native task trace and all six native endpoint statistics.
Its future report is
`benchmarks/work/qfo_corrected_fastoma_assessment_admission_20260918.json`.
Successful process execution alone does not admit an accuracy result.

The comparison exporter now accepts admitted FastOMA reports and preserves
supplied-tree semantics. GO/EC similarity and FAS are not F1. The arithmetic
mean of six endpoint statistics remains a project-defined secondary summary.
No existing comparison table has been overwritten or augmented with pending
results.

## Verification And Limits

The complete unit suite passed: 3,541 passed, one skipped. All 86 focused
scoring/admission/export tests passed, including the added mapping-drift
rejection check. Both new batch scripts pass Bash syntax checking, and the
frozen scoring CLI loads successfully.

An earlier full-suite run caught two unfinished integration cases: the old
generic FastOMA test fixture and the missing export adapter. Both are resolved;
the final suite is green. Synthetic orchestration tests mock native metric
validation; they do not substitute for the queued real-data admission.

At submission, both jobs are dependency-pending. No corrected FastOMA score,
paired uncertainty result, dedicated timing measurement or publication-ready
claim is established by this integration. Shared-host scoring time is not a
matched-resource inference timing result.
