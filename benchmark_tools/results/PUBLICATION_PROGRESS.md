# Publication Progress

Objective: complete the seven-part publication goal, with QfO and OrthoBench
primary and Three Kingdoms supplementary. This ledger is not a claim of
publication readiness. Existing benchmark outcomes are development-exposed.

## Status (2026-09-16)

| Requirement | Status | Evidence / next action |
| --- | --- | --- |
| Baseline and scoring audit | In progress | Historical audit: SCORING_AUDIT_20260910.md; refresh OrthoMCL and correct output semantics before freezing tables. |
| Independent generalization | Not complete | Inventory all development exposure; choose a genuinely independent reference and freeze protocol before scoring. |
| HMM and phylogeny ablations | Not complete | Existing experiments are exploratory; design matched controls and preserve negative results. |
| Uncertainty and error analysis | In progress | Primary OrthoBench comparisons now have paired RefOG bootstrap estimates; other comparisons, QfO uncertainty, error strata, and tracing remain. |
| Robustness and efficiency | Not complete | Validate simulator, multi-seed conditions, tree perturbations, and matched resource measurements. |
| Biological usefulness | Not complete | Prespecify independent families and evidence before selecting examples. |
| Publication/reproducibility package | Not complete | Generate artifacts from audited records; draft manuscript and claim-to-evidence checklist; release/archive remains pending. |

## Verified Current Evidence

- Slurm jobs 20909 (OrthoMCL inference) and 20910 (QfO scoring) both report
  COMPLETED, exit 0:0; elapsed 3-15:16:11 and 00:50:55 respectively.
- The existing QfO OrthoMCL score used cross-species `all_ortho.mtx` edges.
  The installed 1.4 README identifies `all_orthomcl.out` as the final result
  and `tmp/all_ortho.mtx` as the weight matrix. The run's `executeMCL`
  implementation passes that matrix into MCL. Therefore the existing score
  is a pre-clustering diagnostic, not a score of the final clustering.
- The earlier assertion that final-group expansion was necessarily an invalid
  OrthoMCL conversion is withdrawn. Preserve both representations; evaluate
  final group-derived pairs separately, without choosing by observed score.
- Related official documentation also separates potential pairs and clustered
  groups: https://github.com/stajichlab/OrthoMCL/blob/master/doc/OrthoMCLEngine/Main/UserGuide.txt
  (version 2; supporting context, not a replacement for the installed 1.4 code).
- Existing sample-output changes and untracked experimental results predate
  this publication work and must not be reverted or included accidentally.

## Open Audit Questions

- Quantify BLAST query setup failures, affected input proteins, and final
  membership. Then assess reference-family coverage and whether a repair is
  scientifically warranted; a zero process exit code does not resolve this.
- Confirm final versus intermediate output semantics for every retained tool,
  especially OrthoFinder hierarchical versus legacy groups and ProteinOrtho
  graph stage. Do not infer semantics from method labels alone.
- Consolidate commands, checksums, measured resources, historical resumptions,
  and source revisions. Do not overwrite historical provenance with current HEAD.
- Do not call the historical development/validation RefOG split independent
  after its outcomes have been repeatedly inspected.

## Completed Milestone: BLAST Diagnostics

`audit_orthomcl_blast.py` and its regression tests reproduce the audit in
`orthomcl_blast_audit_20260916.json`, with per-input checksums and per-gene records.
There are 53 failed queries (46 statistics failures, seven short queries),
0.0054275% of 976,504 inputs; all 53 map to QfO identifiers and none occur in
the 79,696 final groups. There are also 150 distinct proteins with U-to-X
warnings, of which 138 occur in final groups. The grouped total is 774,272.
Reference-family impact and a justified repair decision remain unresolved.

The separate `run_orthomcl_final_groups_qfo.slurm` workflow preserves the
existing pre-clustering result and refuses to overwrite existing artifacts.
Final-group scoring is required by output semantics, not selected by score.
Job 20916 launched this workflow from commit `bef05f6` on 2026-09-16, with
8 allocated CPUs and 150 GB memory (conversion and scoring only). It is
running all six QfO assessments. Logs:
`qfo_benchmark/scoring/orthomcl_final_groups_20916.log`.
Existing pre-clustering scores are recorded in
`orthomcl_preclustering_qfo_20260916.json` (secondary project mean
0.7244137191145797); do not relabel these as final-group results.

## Completed Milestone: Primary OrthoBench Uncertainty

The protocol is `ORTHOBENCH_UNCERTAINTY_PROTOCOL_20260916.md`; results and
generated report are `orthobench_paired_uncertainty_20260916.json` and
`ORTHOBENCH_UNCERTAINTY_20260916.md`. The retained point estimates reproduce.
Satellite_v2 minus full OrthoFinder F1 is +1.3696 percentage points, with
paired percentile 95% CI [-4.5037, 7.9168]. High sensitivity minus full
OrthoFinder is -2.3775 points, CI [-8.1780, 3.3483]. Neither F1 interval
establishes an advantage. Precision is higher and recall lower for both
OrthoHMM configurations; report the tradeoff, not overall superiority.
These intervals do not account for historical selection on this benchmark.

## Release Checks Still Needed

The authorized remote accepted milestone `bef05f6`. GitHub reported 21
dependency alerts (one critical, seven high, eleven moderate, two low) during
push. This is an untriaged remote notification, not a validated assessment
of runtime exposure. Dependency review is required before release; do not
silently update pinned benchmarking environments and change their provenance.
