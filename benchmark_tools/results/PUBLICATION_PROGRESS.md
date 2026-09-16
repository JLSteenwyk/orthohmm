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

## Ablation Audit And Replay Controls

`PUBLICATION_ABLATION_PROTOCOL_20260916.md` specifies the profile-expansion,
candidate-expansion, and reconciliation factorial design, separate refinement
and membership-filter diagnostics, and the still-required matched sequence
search control. Existing explorations do not substitute for these controls.

The replay audit found that `replay_phylogeny.py` did not pass production
satellite_v2 membership constraints into reconciliation. This does not alter
the full production path, which already passes them. The replay now accepts
and validates the merge trace, records its checksum and selected policy, and
rejects silent omission when a production trace exists next to the supplied
candidate checkpoint. An explicit unconstrained option supports the separate
filter ablation. Historical replays remain labeled as originally executed;
do not retroactively claim they reproduce satellite_v2.

All 8,440 historical OrthoBench production trace records validate against
the preserved candidate partition. The replay and phylogeny-pipeline unit
tests pass (25 tests), including argument propagation, malformed traces,
candidate mismatch, and fail-before-output protection. Frozen YGOB source,
inputs, settings, and queued job were not changed. Latest live check:
`20916` RUNNING at 49m13s; `20917` PENDING for resources.

Next ablation gate: freeze executable commands and input manifests and
demonstrate cached replay equivalence to the production baseline before
launching and interpreting factorial cells. Shared-reference-resource review,
YGOB report wiring, and all previously listed publication work remain open.

## YGOB Reference-Resource Review

`YGOB_REFERENCE_RESOURCE_AUDIT_20260916.md` records the source/command audit,
primary-source citations, checked file hashes, known family overlap, and
limits of the transfer claim. The reviewed frozen inference paths use input
proteomes rather than supplied YGOB groups; the launcher's reference access
is a checksum preflight, not an inference label input. YGOB's sequence-plus-
synteny curation is not wholly independent of sequence-based evidence or
Saccharomyces annotation history. The audit supports the bounded novel-taxon
experiment, not unrestricted independence. Full per-pillar historical
resource ancestry and redistribution permission remain unresolved.

This completes the planned source-level resource review for the YGOB scoring
gate, with those limitations retained. It does not waive inference/output
completion and arithmetic checks, and does not cover all historical tools.
No held-out scores were inspected and no inference settings were changed.

## OrthoMCL Failure Reference Impact And Baseline Decision

`audit_orthomcl_reference_impact.py` and its native Darwin query script
measure direct reference exposure of the 53 previously audited failed
queries. `orthomcl_reference_impact_20260916.json` preserves per-protein
annotations, reference and executable hashes, and native-log provenance;
`ORTHOMCL_FAILURE_IMPACT_20260916.md` explains the baseline decision.

None of the 53 occurs in mapped SwissTrees/TreeFam-A cases or VGNC's 23,934
asserted pairs, and none has an EC annotation. Four have experimental GO
annotations. All 53 have FAS annotation entries, but only 46 have nonempty
feature-type dictionaries. These are direct exposure counts, not a bound on
indirect clustering changes or a measured counterfactual score difference.
The native TreeFam-A file is a single pooled case, not a single gene family.

Retain the unchanged standard OrthoMCL 1.4 baseline with explicit failure
disclosure. No full replacement BLAST run is justified solely by this audit;
altered masking or sequence handling would require a separate diagnostic
configuration, not silent replacement. The modified-search counterfactual
remains unmeasured. Seventeen focused tests pass, and native execution
completed with no reported errors/warnings. Initial audit v1 remains on disk;
v2 clarifies case counts and records actual FAS feature content.

Last live job check: final-group QfO scoring `20916` RUNNING at 1h02m23s;
YGOB inference `20917` PENDING for resources. No active jobs were restarted.

## Generated Historical Accuracy Figures

`plot_publication_accuracy.py` generates three figures directly from the
audited comparison: OrthoBench precision/recall plus supplementary BUSCO F1,
all six QfO native endpoint coordinates, and paired OrthoBench differences
with nominal and multiplicity-adjusted intervals. PNG/PDF/SVG versions and
a coordinate/provenance manifest are in `figures_accuracy_20260916`;
captions and reproduction instructions are in `FIGURE_CAPTIONS_20260916.md`.

Six plotting tests pass, including pending-result exclusion, native-axis
validation, invalid-value rejection, and actual interval coordinates.
All three PNGs were visually inspected for clipping, overlap, visible
pending status, and correct scale labels. QfO uncertainty bars are omitted
because the native recorded fields do not have a uniform interpretation.
The pending OrthoMCL result is not replaced by its pre-clustering diagnostic.
These are explicitly work-in-progress figures, not publication readiness.

## Verified Historical Profile/Refinement Controls

`audit_historical_profile_ablation.py` verifies the original replay cache,
all 12 input FASTA hashes, four stage-partition hashes and membership, and
the corrected fresh production metrics/output. All stages cover 251,378
genes without duplicate membership. The final replay partition is exactly
byte-identical to the corrected fresh production partition (`8ee100f...`).
This comparison uses the self-hit-fixed fresh run, not the superseded earlier
fresh run. It is historical endpoint equivalence, not proof of every current
source branch or intermediate normalized-hit value.

The audited scorer recomputes all four stages in
`historical_profile_ablation_audit_20260916.json`; the generated table is
`HISTORICAL_PROFILE_ABLATION_20260916.md`. F1 is 66.279051 for multipass,
69.763388 after cluster refinement, 66.825390 with profile expansion alone,
and 70.358998 with both. Profile expansion's observed F1 increment after
refinement is +0.595610 points; refinement's increment with profiles is
+3.533607 points. These descriptive controls retain HMM initial search and
cannot establish an overall HMM-versus-sequence-search advantage.

Five new validation tests plus three scorer tests pass. The four-cell
historical audit is reusable evidence, but current-source replay validation,
the eight-cell expansion/reconciliation factorial, QfO counterparts, matched
sequence-search control, and controlled per-arm efficiency measurements
remain required. Latest job check: `20916` RUNNING at 1h12m54s; `20917`
PENDING for resources. Frozen YGOB settings and active jobs were unchanged.
Full unit-suite verification at this milestone: 500 passed in 16.62s.

## Pinned-Source Replay Check Workflow

`run_publication_replay_check.py` validates the historical audit's inputs,
cache, and all four partitions before launching a label-blind replay from
the detached `7f3a9e4` source tree. Current production core has no tracked
differences from that pinned revision. Explicit settings are BLOSUM62,
resolution 0.1, Leiden seed 4, one profile pass, and minimum profile species
one; no jackknife or benchmark-scoring argument is enabled.

The workflow writes source/environment/preflight provenance before inference,
records GNU time separately, preserves failures, and compares all four
resulting partitions with the verified historical partitions. Byte equality
is recorded separately from membership equality so ordering changes cannot
be mistaken for biological differences. Six focused tests pass across the
workflow comparison and historical partition-validation helpers.

Submit with 32 CPUs, an exclusive allocation, and an afterok dependency on
YGOB job `20917`. This is an incremental cached replay, not a controlled
end-to-end runtime. Its successful completion is a prerequisite to reusing
these stages in the new factorial; no new factorial accuracy result has
been evaluated by this workflow.

After committing and pushing launcher `4ae0a82`, a detached launcher worktree
was created at `benchmarks/work/publication_launchers_4ae0a82`. Slurm job
`20919` was submitted with 32 CPUs, 64 GB, a two-hour limit, exclusive
allocation, and `afterok:20917`; `scontrol show job 20919` confirms that
dependency is unfulfilled and the job is pending. It reads the committed
historical audit from the pinned launcher tree and executes inference from
the separate existing `publication_method_7f3a9e4` worktree. Logs will be
`benchmarks/work/ob_replay_check_20919.log`; fresh outputs will be under
`benchmarks/results/publication_ob_replay_check_v1`.

OrthoMCL final-group scoring `20916` remained running at 1h18m13s, and YGOB
`20917` remained queued. The unrelated `foxy_unique720` job was observed but
not modified. Submission is not completion; inspect replay verification
status before using its stage outputs as current-source evidence.

## Completed OrthoMCL Final-Group QfO Assessment

Slurm accounting confirms job `20916` COMPLETED, exit `0:0`, elapsed
01:20:11. All six native endpoint tasks and consolidation completed with
exit zero in `qfo_benchmark/scoring/orthomcl_1_4_final_groups/stats/trace_2026-09-16_10-02-58.txt`.
The final-group workflow metadata records exit zero and a finish timestamp.
The original final-group inputs and converter pass `inputs.sha256` checks;
both generated pair files pass their recorded SHA-256 checks.

The comparison builder now verifies successful workflow metadata, exact
manifest membership, and actual pair-file hashes rather than accepting the
mere existence of a completion manifest. Six added tests cover pending,
success, mutation, failure, unfinished, duplicate, and incomplete cases
(the mutation check shares the success test). All 21 focused comparison,
plotting, and QfO-summary tests pass.

New snapshot `publication_comparison_orthomcl_complete_20260916.json` and
`PUBLICATION_COMPARISON_ORTHOMCL_COMPLETE_20260916.md` preserve the earlier
pending snapshot. OrthoMCL final-group endpoints are VGNC F 0.640957,
SwissTrees F 0.758358, TreeFam-A F 0.717292, EC 0.922140, GO 0.463019,
and FAS 0.729298. The project-defined secondary mean is 0.705177, not the
pre-MCL diagnostic mean 0.724414. No inference or parameter selection was
repeated to obtain this corrected final-output assessment.

The separate scoring GNU-time log records 1:19:48 wall time, 14,762.84 user
CPU seconds, 422.43 system seconds, and 34,127,596 KiB maximum RSS. These
are scoring measurements, not OrthoMCL inference costs or a verified
simultaneous process-tree memory peak. The Nextflow FAS task dominated at
1h16m55s and reports 34.2 GB peak RSS under its own accounting convention.

Regenerated PNG/PDF/SVG figures and provenance manifest are in
`figures_accuracy_orthomcl_complete_20260916/`; all three PNGs were visually
checked for clipping, overlap, and correct labels. The six-panel QfO figure
now includes final-group OrthoMCL. Existing endpoint definitions and
uncertainty caveats are unchanged; these remain work-in-progress figures.

YGOB `20917` remains pending for resources; pinned replay `20919` depends
on its success. No YGOB accuracy was inspected. Concurrent non-Slurm
STRUCTURAL_GENOMICS IQ-TREE processes were observed, so exclusive Slurm
allocation alone cannot establish controlled timing conditions. Unrelated
workloads were not changed. Independent validation, current-source replay,
prospective ablations, robustness, application, and release remain open.

## Frozen YGOB Report Assembly

`report_ygob_validation.py` assembles the four-method score table and the
prespecified paired differences using the independently tested co-membership
scorer. The method set is exact, not inferred from available successful
outputs. Bootstrap settings are fixed at 20,000 PCG64 replicates with seed
20260917; the sequence-only diagnostic is excluded from the two comparisons
and six-metric multiplicity family. Reports label within-species-inclusive
homolog-group recovery, projection of non-reference inputs, coverage, and
the limitations of novel-taxon transfer explicitly.

Its OrthoFinder MCL checkpoint adapter uses the existing native parser and
requires a one-to-one original-ID map matching all inference genes, unique
cluster membership, and complete checkpoint coverage. Unlike checkpoint
validation, final-method scoring permits genuinely missing predictions and
reports coverage rather than silently imputing them.

Eight new synthetic tests cover fixed reporting, the diagnostic exclusion,
method-set validation, and valid/invalid checkpoint mappings. The 21 focused
report/scorer/converter tests pass; the full unit suite passes 515 tests in
18.46s. No held-out accuracy has been computed or inspected. The module is
report assembly, not an end-to-end completion verifier: its output explicitly
marks completion gates as unverified by this module. A caller must still
verify terminal job status, tool completion/versions, source and input hashes,
the frozen reference, and overlap/resource audits before real evaluation.

Latest live scheduler check: YGOB `20917` PENDING Resources, replay `20919`
PENDING Dependency. The frozen inference job, parameters, and input files
were not changed. The next validation milestone is the gated evaluation
entrypoint, followed by real scoring only after successful inference.
