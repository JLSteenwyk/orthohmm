# Publication Progress

Objective: complete the seven-part publication goal, with QfO and OrthoBench
primary and Three Kingdoms supplementary. This ledger is not a claim of
publication readiness. Existing benchmark outcomes are development-exposed.

## Status (2026-09-16)

| Requirement | Status | Evidence / next action |
| --- | --- | --- |
| Baseline and scoring audit | In progress | Historical audit: SCORING_AUDIT_20260910.md; refresh OrthoMCL and correct output semantics before freezing tables. |
| Independent generalization | In progress | YGOB v7 acquired and audited before scoring; taxon/family overlap, reference semantics, and committed evaluation freeze remain required. |
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
  retaining the resolved OrthoFinder and ProteinOrtho checks below. Do not
  infer semantics from method labels alone.
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

## Completed Milestone: Consolidated Comparison

`publication_comparison.py` generates `publication_comparison_20260916.json`
and `PUBLICATION_COMPARISON_20260916.md` from the retained OrthoBench and
Three Kingdoms audits, the primary paired analysis, and the official QfO
assessment JSONs. It preserves QfO axes, participant identity, standard-error
fields, and per-metric checksums, rather than retaining only the custom mean.
It validates Three Kingdoms counts and reproduces primary OrthoBench scores.
Missing QfO results remain pending; the pre-MCL OrthoMCL graph is a separate
diagnostic. Rerun the generator after job 20916 completes to incorporate its
final-group assessment. The report explicitly does not freeze the baseline:
raw-output provenance for some OrthoBench competitors, matched resources,
and the other publication requirements remain open.

ProteinOrtho graph-stage semantics are now confirmed. The installed
`proteinortho_6.3.6--h2b77389_0.sif` contains `/usr/local/bin/proteinortho6.pl`,
which describes `.proteinortho-graph` as the clustered graph (line 383) and
generates it after removing cut edges (line 1852). The
[versioned 6.3.6 manual](https://gitlab.com/paulklemm_PHD/proteinortho/-/raw/v6.3.6/README.md)
agrees under "Clustering Output (step 3)". The retained native-pair choice
therefore does include clustering, unlike the OrthoMCL matrix diagnostic.

## Completed Milestone: Candidate Generalization Data Audit

`INDEPENDENT_VALIDATION_CANDIDATES_20260916.md` documents the current exposure
inventory and candidate-selection rationale. `audit_ygob_overlap.py` audits
the downloaded YGOB v7 snapshot against QfO, OrthoBench, Three Kingdoms, and
test samples. Results and checksums are in `ygob_overlap_20260916.json`.
The audit found 107,277 ON proteins, substantial S. cerevisiae sequence
overlap, and two genes with ambiguous pillar membership. No candidate scores
were calculated. Exact-match absence is not a family-disjointness test.
Dataset acquisition succeeded via the official HTTP endpoint; raw sequences
remain outside Git. Next: settle exclusions and output-level semantics,
complete taxonomic/homology-overlap screening, and commit the evaluation
freeze before launching independent inference.

## YGOB Input Preparation And Scientific Freeze

`prepare_ygob_validation.py` verifies the acquired snapshot and prepares
83,404 proteins from 16 non-Saccharomyces species. The reference has 83,391
genes in 10,250 pillars; 13 genes from ambiguous pillars remain as inference
inputs but are not scored. Input hashes and exclusion lists are recorded in
`ygob_validation_inputs_20260916.json`. Preparation and overlap tests pass.

`YGOB_VALIDATION_PROTOCOL_20260916.md` freezes the scientific specification
before method runs: curated-group recovery, specified micro statistic,
paired pillar bootstrap, two OrthoHMM-versus-full-OrthoFinder contrasts, and
an OrthoFinder MCL checkpoint diagnostic. It explicitly does not claim
family-disjointness. Homology screening and shared-resource checks are still
required before accuracy interpretation. `run_ygob_validation.slurm` performs
only inference from pinned production code, not scoring.

The scientific freeze was committed and pushed as `ec92413` before submitting
Slurm job `20917`. It requests an exclusive allocation, 32 threads per method,
128 GB memory, and sequential fresh high-sensitivity, satellite_v2, and full
OrthoFinder runs. Frozen OrthoHMM source is the detached worktree at `7f3a9e4`.
The job is queued; do not interpret pending status as failure or restart it.
Log path: `benchmarks/work/ygob_validation_v1/inference_20917.log`.
Complete the prespecified homology screen while waiting and before scoring.

## YGOB Overlap Gate And Scoring Arithmetic

The homology-screen implementation was committed and pushed as `fc42cf1`.
Slurm job `20918` runs the frozen DIAMOND screen with 32 CPUs against QfO,
OrthoBench, and Three Kingdoms development inputs. Job `20917` now has an
`afterok:20918` dependency, enforcing successful screen completion before
validation inference. At the latest check, `20918` was running its search,
`20917` was pending on that dependency, and final-group OrthoMCL QfO scoring
job `20916` was still running. No replacement jobs were submitted.

`score_ygob_groups.py` implements the frozen group co-membership statistic,
strict native-group readers, explicit input-ID validation, scoring-universe
projection, exact recovery, and coverage counts with stated denominators.
Singletons and within-species pairs are included as specified. Cross-pillar
false positives are allocated half to each incident pillar for the paired
bootstrap; batched sampling avoids a full 20,000 by 10,250 count matrix.
Ten tests check explicit pair enumeration over 100 random partitions,
excluded genes, missing genes, malformed membership, adapters, paired
resampling arithmetic, and batch invariance. Together with the screen and
input-preparation tests, 20 tests pass. No YGOB method scores were inspected.

Next: complete the gated runs and shared-reference-resource audit, then
connect verified output manifests to the scorer and generate the frozen
two-contrast validation report. The OrthoFinder sequence checkpoint remains
diagnostic and must not enter the six primary/secondary contrast-metric
multiplicity count. Publication ablations, simulations, resource comparisons,
biological application, and manuscript deliverables remain open.

The subsequent scheduler check confirms `20918` COMPLETED with exit `0:0`
in 5m23s. The saved report `ygob_homology_screen_20260916.json` contains
input/tool/source provenance and hit-file checksums, verified against disk.
71,714/83,404 proteins (85.98%) and 6,952/10,250 pillars (67.82%) have a
qualifying development hit. This establishes substantial family overlap;
novel-taxon transfer remains the intended claim, not family-disjointness.
No families or endpoint definitions were changed in response to the screen.
`20917` now waits for resources after its dependency succeeded; `20916`
continues running. Full unit-suite verification: 469 passed in 18.73s.
