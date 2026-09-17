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

## Label-Blind Inference Verification Entry Point

`verify_ygob_validation.py` queries Slurm accounting for the exact parent
job, requiring COMPLETED and exit 0:0 before opening outputs. It checks
runner completion metadata, the frozen source revision, tracked source
cleanliness, launcher and OrthoFinder-entrypoint checksums, frozen input
manifest agreement, recorded input/reference hashes, both OrthoHMM harness
completion records and source/input/output manifests, required native output
presence, identical OrthoFinder FASTA copies, and its GNU-time exit status.
Manifest paths must remain inside their expected base directories, without
duplicates; both file lengths and hashes are checked.

This is an inference file/completion verifier, not a claim that every
scientific scoring gate has passed. Its output explicitly leaves exact
command/version review, overlap/reference-resource audit verification,
native ID conversion, and independent reference reconstruction/scoring
open. Successful full-run integration remains untested until completed
inference is available. No accuracy is computed by this entrypoint.

Thirteen new tests cover scheduler ambiguity/failure, early rejection before
output reads, artifact mutation, bytes/hashes, duplicate/escaping paths, and
metadata parsing. All 31 focused verifier/report/scorer tests pass. Executing
the real CLI for pending job `20917` returned exit 1 at the scheduler gate,
before opening prediction files or writing verification output. This is an
expected refusal, not a failed inference or grounds to restart the job.
The prior report assembly milestone remains intact; no held-out outcomes
were inspected and no queued inference configuration was changed.

## Initial Manuscript And Claim Audit

`PUBLICATION_MANUSCRIPT_DRAFT_20260916.md` now provides the study objective,
benchmark/output semantics, statistical methods, prospective validation and
ablation descriptions, evidence-backed development Results, limitations,
and an explicit unfinished data/code availability statement. It does not
claim held-out success, overall superiority, controlled speedups, or a
completed archive. The historical source records are distinguished from
the prospective source pin. Literature bibliography and detailed final
algorithm description remain to be completed alongside the open experiments.

`PUBLICATION_CLAIMS_20260916.md` maps proposed claims to evidence and audits
all seven original work packages. A linked protocol is explicitly not
treated as a completed experiment. Missing provenance, matched HMM controls,
error strata/tracing, simulations, robustness/scaling, biological application,
release/security work, and archival steps remain visible completion gates.

Validated all 29 local evidence links across the two documents and checked
nine principal manuscript point estimates directly against the committed
comparison JSON. Remaining detailed numbers were transcribed from their
linked audited result reports; this draft is not a replacement scoring
pipeline. No inference or scoring code changed at this milestone.
Jobs `20917` and `20919` remain queued as resource/dependency waits; unrelated
workloads and the frozen validation configuration were not modified.

## Evolutionary Simulator Compatibility Test

Inspected the primary ALF source/manual and publication, pinned the author
repository at `b674ab10018c3c0fcc0806434dddf7b36136e2c1`, and ran a small
protein-evolution smoke-test driver against the existing Darwin container.
The unmodified source/engine combination failed before parameter loading
with `Bad VectorABQMode`, exit 1. No simulation or accuracy result exists.
Driver, test parameters, native failure log, source/container hashes, and
reproduction instructions are retained in `SIMULATOR_PREFLIGHT_20260916.md`.
This resolves an infrastructure feasibility question but does not fulfill
the simulation requirement. A compatible official ALF engine or another
published simulator must be evaluated next; no algorithm patches were made.

An overly broad container-environment diagnostic exposed credentials in
tool output. No environment dump or credential values were copied into
repository files, reports, or commits. The user was notified to rotate the
affected credentials. Subsequent diagnostics used targeted queries and the
simulation invocation used `--cleanenv`. This security issue must not be
confused with the unrelated pre-existing dependency alerts.

## Reproducible Zombi Execution Path

The alternative Zombi source is pinned at
`8db13ee4ba007f46c17f38586d31e5aa617c1647`. An isolated overlay environment
adds ETE3 3.1.3 and uses Pyvolve 1.1.0 without changing the base environment.
Source inspection found that Zombi's global seed does not seed Pyvolve's
per-call generator. `run_zombi_seeded.py` therefore supplies a stable
128-bit SHA-256-derived family seed through the native Pyvolve seed API,
preserving other arguments and explicitly supplied seeds. The adapter and
its version are part of simulation provenance, not hidden upstream changes.

`smoke_zombi.py` ran T/G/S for two equal seeds and one independent seed.
All nine stages completed; all 84 products match byte-for-byte between
equal-seed runs, and the independent seed changes sequences. The two
distinct histories include one versus two duplication events; the second
also includes a loss. All outputs, parameters/defaults, commands, source
pin and package versions are recorded in `zombi_smoke_20260916.json` and
explained in `ZOMBI_SMOKE_20260916.md`. Raw v1 (32-bit seed adapter) and
v2 (retained 128-bit adapter) outputs are separately preserved.

Five focused adapter tests pass; the full unit suite now passes 533 tests
in 22.51s. This establishes reproducible simulator execution, not validated
orthology truth or method accuracy. Extant-only sequence extraction,
event/tree/XML truth checks, scientific multi-seed conditions, indel/fragment
limitations, matched method runs, and robustness/resource analyses remain
open. No held-out YGOB outcomes or method parameters were inspected/changed.

## Event-Derived Simulation Truth

`zombi_truth.py` now cross-checks native event histories against reconciled
XML, extant genomes, pruned gene trees, and protein sequences. It exports
extant-only FASTAs with globally unique IDs and event-derived ortholog pairs,
not family cliques. Unsupported transfer/origination histories, inconsistent
outputs, invalid graph structure, duplicate IDs, and invalid proteins fail
closed. Root-origin homolog families are explicitly not a root-HOG truth set.

Both independent smoke histories pass: seed 20260916 has 41 extant genes
and 63 cross-species ortholog pairs; seed 20260917 has 44 genes and 66 pairs.
Reports and hashes are in `zombi_truth_repeat_a_20260916.json` and
`zombi_truth_independent_20260916.json`; interpretation and limitations are
in `ZOMBI_TRUTH_VALIDATION_20260916.md`. Fifteen new tests, including
duplication timing, co-orthologs, loss, and corrupted-output fixtures, pass;
20 tests pass together with the seed-adapter tests.

No method accuracy has been computed on these inputs. Extinct/single-tip
family cases, scientific conditions and multiple-seed evaluation, missingness
transformations, root-time group semantics, and matched method runs remain
open. YGOB and replay jobs were still queued at the start of this milestone;
their inputs and frozen configurations were not changed.

## Prospective Simulation Panel Specification

`PUBLICATION_SIMULATION_PROTOCOL_20260916.md` freezes a ten-seed, seven-condition
design before method inference or scoring: baseline, divergent, turnover,
divergent_turnover, 20% gene missingness, clade-thinned taxa, and a taxon-count
control. Forty native simulations plus thirty derived datasets give seventy
condition/seed evaluations. Genome-level rates, realized-event reporting,
native-history reuse checks, deterministic label-independent transformations,
matched four-CPU method settings, event-derived pair endpoints, seed-level
uncertainty, multiplicity and failure handling are specified explicitly.

This is a scientific protocol, not an executed or fully materialized panel.
An executable manifest with expanded defaults, hashes and exact commands,
tested transformation/scoring/aggregation code, and frozen dependencies is
required before launch. Tree-error/parameter-neighborhood experiments,
representative scaling, curated validation and biological application remain
separate requirements; the small synthetic panel does not replace them.

The native source confirms that fully lost single-node trees may have no
sequence file. The truth adapter now permits that only when independent
event/XML/genome/pruned-tree checks establish zero survivors. Missing sequence
files for survivors still fail. Two new fixtures cover completely lost and
single-survivor families; 17 truth tests pass. Re-evaluating both real smoke
histories confirms unchanged family membership and ortholog pairs (comparison
normalizes Python tuples to JSON lists). The historical reports retain their
original source provenance and were not overwritten.

## Implemented Simulation Transforms And Pair Scoring

`simulation_conditions.py` implements the frozen missing20, uneven-clade,
and taxon-count control selections, using no labels or sequence information.
It projects event-derived truth and scores all retained-universe predicted
pairs, including cross-family false positives. Invalid IDs/pairs and
duplicate truth are rejected; orientation/duplicate prediction rows are
explicitly canonicalized. Pair-endpoint coverage and undefined ratios are
labeled. Successful method completion remains a separate caller gate.

`derive_simulation_conditions.py` validates native baseline truth before
exporting transformed FASTAs, projected truth and provenance. Existing
outputs are refused, and inapplicable tree selections remain inapplicable.
The real four-species smoke export retained 33 genes/38 true pairs for
missing20, 31/32 for uneven_taxa, and 30/30 for taxon_count_control, from
41/63 baseline. This is a workflow check, not a scientific panel outcome.
The manifest is `zombi_transforms_smoke_20260916.json`; interpretation and
reproduction are in `SIMULATION_TRANSFORMS_20260916.md`.

Thirteen new tests cover selections, false-positive accounting and exports;
30 focused transform/export/truth tests pass. No method accuracy was
computed, no validation outcomes were inspected, and no frozen method
configuration was changed. Executable panel materialization, seed-level
aggregation and matched method execution remain the next simulation steps.

## Simulation Manifest Builder

`prepare_simulation_panel.py` expands all pinned native defaults into exact
T/G/S parameter files for the forty prespecified configurations and indexes
all seventy native/derived datasets. It records workflow/native-source hashes,
default and generated parameter hashes, exact generation/truth/transform
commands, the ten seeds, matched-history checks, six resolved simulation
dependency versions, and a package-version inventory without direct URLs or
environment-variable dumps. Existing panel directories are refused.

Four focused tests verify dimensions, matching T/G parameters within the
divergence contrasts, all overrides, default preservation, and rejection of
out-of-protocol seed/condition requests. The builder executes no simulator
or method inference. Pin it in a detached worktree before materialization;
generation-runner integrity checks, seed aggregation and the method launch
manifest remain open gates. A version inventory is not yet proof of a
portable environment rebuild or a complete release dependency lock.

Builder `4c5c7b0` was committed/pushed and checked out in detached worktree
`benchmarks/work/publication_simulation_4c5c7b0`. Running it with the isolated
simulation Python materialized `benchmarks/work/publication_simulation_panel_v1`.
The committed manifest is `publication_simulation_manifest_20260916.json`;
`simulation_environment_versions_20260916.txt` preserves the version inventory.
Commands reference the pinned worktree rather than the mutable main checkout.

Verified all 120 generated parameter file lengths and hashes, all recorded
workflow hashes, 70 unique condition/seed dataset entries, and an empty native
output directory. The manifest contains 40 native configurations, 170 staged
generation/truth/derivation commands, and 20 required history-equivalence
checks. Status is `materialized_not_executed`; no successful simulation or
method outcome is inferred from these files. Generation runner checks,
seed-level aggregation and pinned method conversion/launch remain open.

## Prespecified Simulation Seed Aggregation

`summarize_simulation_panel.py` implements the frozen ten-seed/seven-condition
analysis for the two OrthoHMM modes, full OrthoFinder and its diagnostic
checkpoint. It requires all 280 explicit terminal method/dataset rows;
missing or pending rows are not silently considered failures. Successful
scores require consistent counts/metrics, explicit undefined-ratio flags,
and a shared truth hash/input universe among methods for each condition/seed.
Failures/inapplicable rows require reasons and cannot carry imputed scores.

The statistic is the mean of seed-level F1/P/R, not pooled gene-pair counts.
Paired complete-seed resampling uses 20,000 PCG64 draws and seed 20261031,
reset per contrast so identical included seed sets share draws. F1 intervals
carry the prespecified 14-comparison Bonferroni adjustment; P/R intervals are
exploratory. Reports retain included/excluded seed IDs, all original rows,
available-case means, failure fractions and conditional-estimation caveats.
Zero/one complete pair sets do not produce misleading bootstrap intervals.
CLI output records source/input hashes, command, Python and NumPy versions.

Eleven new tests verify direct multinomial resampling, the distinction from
pooled counts, diagnostic exclusion, failures, invalid provenance/metrics,
missing rows and insufficient seeds. Twenty-six focused aggregation/scoring/
parameter tests pass. Only synthetic test records were evaluated; no real
method scores or held-out validation outcomes were inspected. The generation
runner integrity checks and pinned method execution/conversion manifest
remain the outstanding panel launch prerequisites.

## Hash-Gated Simulation Generation Runner

`run_simulation_generation.py` executes one named run from the immutable
panel manifest. Before execution it verifies the externally supplied
manifest SHA-256, exact simulator revision, source/default/workflow hashes,
all 120 parameter hashes, interpreter version and complete recorded package
inventory. It rejects escaping artifact paths, changed inputs, ambiguous
run labels, and any existing stage output; there is no automatic restart
that could erase expensive or failed evidence.

Each stage has a native log, separate GNU-time log, explicit argv, start/end
times and exit status. The first failed stage stops later stages while
preserving logs, partial native output and manifest/runner provenance.
Successful stage exits and successful output-hash inventory recording are
separate fields. Inventory failures are retained as terminal failures, not
silently left as verified completion. The runner never executes methods or
computes accuracy; native truth validation is an explicit manifest stage.

Four new tests cover artifact corruption, path escapes, manifest mismatch,
failure preservation/stopping, non-overwriting and successful stage records.
Nineteen focused runner/builder/aggregation tests pass. The real check-only
invocation passed for `baseline_20261001`, against manifest SHA-256
`ee31ea38d3b5c047abf80f04636f06838959f941d6704649a216bb93958343b2`.
No scientific-panel simulation was launched. The next gate is the pinned
method execution/conversion manifest and remaining matched-history checks;
all existing validation and replay jobs were preserved.

## Native Simulation Method Adapters

`simulation_method_outputs.py` connects high-sensitivity groups, native
OrthoHMM phylogenetic pairs, full OrthoFinder tables and its sequence-only
MCL checkpoint to the simulation pair scorer. Native species annotations
must agree with the input mapping. Artifact discovery requires unique paths.
OrthoFinder default output must contain every directed species-pair table,
including empty tables, and opposite orientations must agree. It then uses
the existing audited native converter. The MCL adapter reuses the strict
complete-universe ID restoration tested for YGOB. Unknown/duplicate IDs,
missing tables and inconsistent outputs cannot become empty successful
predictions. Tool completion remains a separate execution gate.

Nine new format/completeness tests pass; 31 focused adapter/converter/scorer
tests pass. A real unscored OrthoFinder 3.1.5 run used the four-species,
41-gene simulation smoke FASTAs with `-t 4 -a 4 -S diamond`. Native execution
and GNU time record exit zero; elapsed wall time was 5.80 seconds under
uncontrolled machine load, not matched efficiency evidence. Full output
passed all 12 directed table checks; checkpoint restoration also passed.
Both conversions yielded 63 pair rows, without computing overlap/F1 against
truth. Equal pair counts alone are not proof of equal predictions or accuracy.

Input/native artifact hashes, converter hashes and execution-log hashes are
in `orthofinder_simulation_adapter_smoke_20260916.json`. Raw inputs/results
and logs are preserved at `benchmarks/work/orthofinder_simulation_adapter_smoke_v1`.
This smoke dataset is not one of the 70 scientific-panel entries; no panel
or held-out YGOB outcome was inspected. Pinned method execution commands,
completion checks and matched-history validation remain to finish before
scientific panel launch.

## Frozen Simulation Method Command Builder

`prepare_simulation_methods.py` constructs exact inference commands for all
70 prespecified datasets from the immutable generation manifest. OrthoHMM
uses source `7f3a9e40dd7e79f842cc2c11fb8b548f9a802806`, explicit high-sensitivity
settings and four CPUs/worker threads; satellite_v2 has explicit inferred
tree/rooting/pair settings. OrthoFinder is checked through installed package
metadata as version 3.1.5 and uses four search/analysis threads plus DIAMOND.
Its sequence-only checkpoint has a parent relation, not an extra command
or independent timing. OrthoFinder receives a separate, verified input copy.

The builder records tracked OrthoHMM sources, adapter/scorer sources, resolved
tool entrypoints, the installed OrthoFinder distribution files, both Python
package inventories and required environment overrides. Reference truth paths
are kept out of inference argv. All external-tool auxiliary dependencies and
actual runtime resolution still need explicit accounting; the manifest does
not claim complete portable reproducibility from entrypoint hashes alone.

Two new command/gate tests and nine native-adapter tests pass. Installed
OrthoFinder metadata confirms 3.1.5 with 189 distribution entries. Pin this
builder in a detached worktree before producing the committed command
manifest. No panel inference or held-out accuracy has been evaluated.

Builder `213b180` was committed/pushed and pinned in detached worktree
`benchmarks/work/publication_methods_frozen`. It generated
`publication_simulation_methods_20260916.json` with SHA-256
`4f0717fa8b0c34d5a6982a3193772d64bc39a3d75b8fe23195d9ce1b21eee186`.
Verified all 186 recorded source/tool artifacts against actual lengths and
hashes, 210 inference commands, and 70 checkpoint definitions. Command and
adapter references point to pinned worktrees; no mutable main source is
used for the prospective OrthoHMM core or adapters.

Resolved tool paths are the existing MAFFT 7.525-with-extensions launcher,
FastTree_v220 binary, diamond-linux64 binary, and OrthoFinder 3.1.5 virtual
environment launcher. Exact paths/hashes are in the manifest. Its output
root `benchmarks/results/publication_simulation_methods_v1` does not exist:
no scientific-panel inference has started. Remaining execution work includes
matched-history verification, generation scheduling/completion, and a method
runner that verifies generated input hashes and terminal success before
native conversion/scoring. These are execution tasks, not an invitation to
change the frozen scientific settings after seeing results.

### Simulation generation scheduling (2026-09-16)

Added `run_simulation_generation.slurm`: 40 native generation tasks, capped
at two concurrent tasks, each with one CPU, 8 GiB requested memory, and a
12-hour limit on bizon. Manifest and pinned runner hashes are checked before
task selection; the existing runner checks parameters, sources, environment,
and absent output paths. It preserves failed artifacts instead of restarting.
The frozen settings and 70-dataset design are unchanged. Generation is not
method inference, and these shared-machine timings are not scaling evidence.

Shell syntax and the last array entry's full check-only preflight passed;
index 40 was rejected. Eight generation runner/builder tests passed. Commit
and push this launcher before submission. YGOB job 20917 remains pending for
resources; replay job 20919 remains dependent on it; unrelated job 20915 is
running and has not been modified. Next: record submission, validate paired
histories after generation, and execute the frozen method comparisons.

Launcher milestone `7987474` was pushed before submission. Slurm accepted
array job **20920** (indices 0-39, concurrency 2), with logs under
`benchmarks/work/publication_simulation_panel_v1/slurm/20920_%a.log`.
The first eight tasks completed native stages and truth export successfully;
both completed baseline tasks also finished their derived conditions.
Remaining tasks are running/queued, not assumed successful.

Added `verify_simulation_histories.py`, which requires terminal generation
success, exact manifest and stage commands, and verified output checksums
before comparing every T/G biological file. Only the two copied parameter
files are excluded. Missing histories and unrecorded files fail validation;
history mismatches are retained in a report and return failure. Full-panel
verification requires all 40 runs and the 20 prespecified comparisons.
Three new tests plus four generation tests passed. The existing deterministic
smoke pair matches all 61 history files. The first scientific baseline and
divergent pair (seed 20261001) passes completion/provenance checks and matches
all 437 history files. No accuracy outcomes were evaluated. Method execution,
full-panel validation, other robustness experiments, ablations, and the
remaining publication requirements are still outstanding.

Full unit-suite verification after this change: **596 passed in 23.46s**.
Slurm accounting independently confirms array tasks 0-9 completed with exit
`0:0`; task-generation timings are approximately 19-21 seconds, not inference
benchmarks. Other tasks remain in progress.

### Frozen simulation inference executor (2026-09-16)

Added `run_simulation_methods.py` and a 70-task Slurm launcher. These execute
the already frozen 210 method commands, not a revised scientific design.
Preflight verifies method/generation manifest hashes, recorded core/tool/
adapter files, interpreter package inventories, generated file checksums,
exact species input file sets, and the parent's paired biological history.
OrthoFinder receives a verified copy. PATH resolution is recorded explicitly
as lookup evidence, not a trace of every spawned executable.

Each method's command, exit status, GNU time log, output hashes and failures
are preserved separately. A failed method does not suppress the remaining
methods. Successful exits remain `process_succeeded`, with dataset status
`finished_pending_native_validation`: no native validation or accuracy is
implied. No existing output or evidence is automatically overwritten.
Inapplicable derived conditions require verified generation evidence and
run no inference. Full native validation/conversion/scoring is still needed.

Five new tests cover manifest integrity, copies, failed-method continuation,
inapplicable conditions and missing species files; these and three history
tests pass. Real baseline_20261001 preflight passes without inference.
Generation array 20920 has progressed to tasks 28-29 running, with earlier
tasks finished. Pin and push the executor before scheduling method tasks
behind generation completion; requested resources are four CPUs, 16 GiB,
12 hours per task and at most two concurrent tasks. These shared-machine
timings are descriptive, not controlled scaling evidence.

Executor milestone `66dd9c4d62e9495e104d880578084563eeb54ff8` was pushed
and pinned in `benchmarks/work/publication_execution_v1`. Its Slurm launcher
also passed a real check-only run for array index 0 (missing20_20261001).
Submitted method array **20957**, indices 0-69, with `afterok:20920` so no
method task starts before the full generation array succeeds. Logs reside
under `benchmarks/work/publication_simulation_panel_v1/method_slurm/`.
The executor verifies the required paired history itself before each dataset.
Full unit suite: **601 passed in 22.31s**. At this checkpoint generation
tasks 36-37 are running and method tasks remain pending on the dependency.
Native output validation and scientific scoring remain separate next steps.

### Complete generation and native numerical failure audit (2026-09-16)

All 40 generation tasks are confirmed `COMPLETED/0:0` by Slurm accounting.
All 70 dataset truth files exist. Full paired-history verification succeeded:
20 matched comparisons, 437 biological files per side. Committed inventory
report: `simulation_histories_20260916.json`.

Added native output validation with six tests, including rejection of
nonfinite MCL graph weights and native failure despite process exit zero.
Both baseline_20261001 OrthoHMM runs pass completion/provenance checks, but
OrthoFinder's nominally successful process has no native completion marker
and reports no usable species-tree alignments. Its graph contains nan/inf.
This is why process-success records were deliberately not treated as scores.

Reproduced the failure by calling installed, unmodified OrthoFinder 3.1.5
normalization on saved species 0 versus 1 DIAMOND hits. All 98 raw hits have
the same length product (90,000); fitted coefficients cause overflow, and
all 98 normalized values are nonfinite. Actual search hits exist, so the
initial warning alone was not evidence of a search-program failure.
See `SIMULATION_NATIVE_FAILURE_AUDIT_20260916.md` and the hashed diagnostic
JSON/script. This finding changes the next scientific action: preserve this
fixed-length stress panel, and freeze a separate heterogeneous-family-length
panel before evaluating additional outcomes. Do not reinterpret this failure
as ordinary competitor accuracy or as general superiority of OrthoHMM.

Method array 20957 remains active; original outputs and configurations are
preserved. No accuracy scores have been computed. Native admission checks,
proper explicit failure records, and the heterogeneous-length protocol are
next, alongside the queued YGOB/replay and other outstanding goal work.

Verification: **607 unit tests passed in 22.81s**; fifteen focused native
validator/adapter tests passed. The installed-tool numerical diagnostic
completed successfully and preserved its warnings/counts without modifying
native outputs or the competitor installation.

### Prospective variable-length protocol and successful smoke (2026-09-16)

Protocol `PUBLICATION_VARIABLE_LENGTH_PROTOCOL_20260916.md` was committed
and pushed as `cec7da4` before the smoke test. It specifies new scientific
seeds 20261101-20261110 and a fixed SHA-256 mapping of root-family IDs to
100-500 amino-acid lengths; it is not an empirically fitted proteome model.
All other scientific conditions and frozen method settings remain unchanged.
Keep the original constant-length panel as a distinct stress test.

Native Zombi Sf mode switches to a codon model, so it cannot isolate length
variation while retaining WAG. Added an optional length-aware wrapper around
ordinary S mode: verify the mapping against the frozen rule and exact native
family-tree set, set sequence size for one family call, and restore it in a
finally block. Native evolution and name correction remain unchanged. Default
fixed-length execution is unchanged, and already-running jobs use their
previous pinned source. Added eight tests for stable length assignment,
invalid inputs and wrapper behavior; seventeen driver/builder tests pass.

Ran `smoke_zombi_variable_lengths.py` with the prespecified non-panel seed
20260918, eight species and 100 root families, in
`benchmarks/work/zombi_variable_length_smoke_v1`. Both repeats completed all
native stages with byte-identical products; native event/XML/tree/FASTA truth
validation succeeded for 816 extant genes. Exported proteins have 88 distinct
lengths, minimum 102 and maximum 495, each exactly matching its assigned
family length. OrthoFinder 3.1.5 completed using the unchanged four-thread
command; native version/command/completion checks and finite-graph checks
passed. No smoke accuracy was computed. The complete hash/command evidence
is `zombi_variable_length_smoke_20260916.json`.

Next: materialize and pin the additional panel's length mappings, generation
and method manifests; extend existing verification/aggregation interfaces
without changing original-panel semantics; run the new panel; finish native
failure accounting and scoring for the original stress panel. Independent
validation, ablations, other robustness tests and remaining publication work
are still required.

Full unit-suite verification after the optional length adapter: **615 passed
in 24.29s**. This is implementation and simulator smoke evidence, not a
scientific-panel accuracy result or a publication-readiness claim.

### Materializing the variable-length panel (2026-09-16)

The existing generation builder now accepts an explicit `variable_length_v2`
variant with the new protocol/seeds. Default fixed-length settings remain
unchanged. The new variant writes one deterministic 100-family mapping per
seed, records each mapping as a hashed extra input, and passes it to both
sequence generation and truth export. The generation preflight verifies all
extra inputs before any stage starts. The truth exporter checks every extant
sequence against the frozen length rule and records mapping/rule provenance.

Added a generic pinned-worktree generation launcher; existing v1 launchers,
manifests and running jobs are untouched. Two new tests cover disjoint seeds,
unchanged evolutionary settings and rejection of invalid exported lengths.
All 27 focused parameter/generation/truth tests pass. The new export gate also
passed on the real 816-gene variable-length smoke run, using a new output
directory rather than replacing its earlier evidence. Pin and push the
builder/runner before materializing or launching the scientific v2 panel.

Pinned generation builder/runner `08f2888125c46f376d5e0da2181d18ce067f8d5e`
in `benchmarks/work/publication_variable_generation_v2` after pushing it.
Materialized `benchmarks/work/publication_variable_simulation_panel_v2` with
40 native runs, 70 datasets, 120 parameter files and ten hashed family-length
mappings. Committed-manifest candidate
`publication_variable_simulation_manifest_20260916.json` has SHA-256
`806aa1e5f6976c323ff2f6641dd88e25264e7417749e7d77565294666aee776b`.
Real first-run preflight passed without execution, including extra-input
checks. No scientific v2 native run has started at this checkpoint.

Extended method preparation/execution to accept an explicit frozen manifest
hash, retaining the original hashes as defaults for v1 callers. Generation
provenance is derived from the hash-protected method manifest instead of a
v1-only constant. Native validation sources are included in the new method
manifest. A generic pinned method-array launcher uses the same four-CPU,
two-concurrent-task resource settings. Ten focused method/variable-length
tests pass. Pin this code and materialize/freeze method commands before
submitting v2 generation, as required by the prospective protocol.

Pinned method builder/executor `f5f4e1b662118e7624746546e81505ddeabf9bb8`
in `benchmarks/work/publication_variable_methods_v2`. It generated
`publication_variable_methods_20260916.json`, SHA-256
`f68b0caf508fde8f40d311b5d303f7569a118bd15e8cf960d227a7d1cfee5984`.
All 70 datasets have the unchanged three inference commands and diagnostic
OrthoFinder checkpoint definition. Verified exact equality of the frozen
core revision, tool-entrypoint records and scientific command options with
v1; all recorded source/tool files and both interpreter inventories passed
preflight. The generic generation launcher also passed its full check-only
preflight for the final task (divergent_turnover_20261110).

After committing/pushing both manifests, submitted v2 generation array
**21009** (40 tasks, one CPU and 8 GiB each, concurrency two) and method
array **21010** (70 tasks, four CPUs and 16 GiB each, concurrency two), with
`afterok:21009` on the method array. Both use the pinned worktrees/commits
above; neither changes v1 job 20957. Logs are under the v2 panel's `slurm/`
and `method_slurm/` directories. The scheduler initially reports generation
pending on priority; completion is not assumed. Full unit suite: **618 passed
in 23.00s**.

Next: monitor authoritative terminal outcomes; verify all v2 paired histories
and exported length checks; finish native numerical/completion failure
accounting and scoring for both panels; extend the prespecified aggregation
to the new seed/bootstrap constants without pooling panels. The queued YGOB
validation and OrthoBench replay, HMM/phylogeny ablations, other robustness
and efficiency experiments, biological application and final publication
package remain outstanding.

### Terminal result assembly and panel-specific inference (2026-09-16)

Added `assemble_simulation_results.py`. It requires all 70 scheduler tasks
to be uniquely terminal before any panel scoring, binds execution evidence
to raw job/task IDs and pinned executor sources, revalidates generated inputs
and paired histories, and checks native method completion. Frozen conversion
and pair-scoring source hashes are verified separately from the admission
gate code. Resource logs and prediction artifacts retain checksums.

Only explicitly identified native failures (missing completion or nonfinite
weights in verified output) become native-failure records. File integrity,
command, version and input mismatches raise errors instead of being silently
classified as poor method performance. Execution/interruption failures have
distinct reasons; no failed method receives a score. The OrthoFinder MCL
checkpoint inherits the parent's native admission requirement and has no
independent timing. Successful predictions are scored across the complete
input universe, retaining cross-family false positives.

Extended seed aggregation with explicit fixed_length_v1 and variable_length_v2
protocol choices. Each retains its ten seeds, separate results and fixed
bootstrap seed (20261031 or 20261130); mixed/wrong-panel seeds are rejected.
The old default and paired statistic remain unchanged. Eight new tests cover
terminal-state gating, failure-versus-integrity distinctions, source/job/input
binding, cross-family false positives, and the new bootstrap specification.
All 25 focused tests and **626 unit tests (23.84s)** passed.

Live full-panel CLI preflight correctly refused unfinished task 20957_63
without writing results. Label-blind admission checks on completed task
20957_0 verified scheduler/source/input provenance: both OrthoHMM runs were
admitted and both OrthoFinder views were rejected for missing native
completion despite exit zero. No accuracy was computed in this check.
Commit/pin the assembler before full-panel scoring once the jobs finish.

Assembler `a58fd3a` was committed/pushed and pinned in
`benchmarks/work/publication_simulation_scoring_v1`. Once all v1 tasks were
terminal, it produced `simulation_fixed_length_results_20260916.json` with
280 explicit outcomes. The generated Markdown table includes completed,
failed and inapplicable denominators alongside available-case means and
paired intervals. A renderer test proves all-failed means remain NA. Full
suite after the renderer: **627 passed in 22.33s**.

Actual v1 admission: high sensitivity 70 complete; satellite_v2 64 complete
and six execution failures; OrthoFinder full 70 rejected (48 missing native
completion, 22 nonfinite graphs), with all checkpoint views also rejected.
Thus no OrthoFinder comparative F1 or CI can be estimated on this stress
panel. Baseline descriptive OrthoHMM F1 means are 96.29% and 99.68%, each
from ten successful seeds. Do not interpret these as competitor superiority.
All six satellite_v2 native logs identify insufficient connected single-copy
family coverage for species-tree inference on divergent seeds 20261003,
20261006 and 20261007. Its seven-seed divergent means are conditional and
not directly paired against the ten-seed high-sensitivity means.
See `SIMULATION_FIXED_LENGTH_INTERPRETATION_20260916.md` for scope and paths.

V2 generation 21009 finished, and dependent method array 21010 is running.
Pinned history verification succeeded for all 20 pairs; inventories are in
`variable_simulation_histories_20260916.json`. Re-read all 40 native prepared
exports with Bio.SeqIO and re-applied `check_family_lengths` to 33,618 sequence
instances: every length matched its frozen assignment and exported validation
metadata, and assignment hashes matched. No v2 accuracy has been evaluated.
The broader publication goal remains active; queued YGOB/replay, ablations,
additional robustness/efficiency/application work and final packaging remain.

### Unblocking the OrthoBench replay (2026-09-16)

Audited job 20919's actual batch script and pinned replay implementation:
it reads no YGOB output and does no accuracy scoring. The original dependency
on 20917 and exclusive allocation were scheduling precautions, not required
scientific inputs. This cached correctness check is already excluded from
controlled end-to-end runtime comparisons. Added an explicit scheduling
amendment to the ablation protocol, preserving CPU=32, memory=64 GiB, two-hour
limit, exact command, source revisions, input/cache hashes and output path.

Attempted `scontrol update JobId=20919 OverSubscribe=YES Dependency=`;
Slurm refused permission and the job remained unchanged/pending. Cancelled
only this unstarted job; accounting confirms `CANCELLED by 1000`, runtime
00:00:00. No replay output directory exists. Recovered the original submitted
script with `scontrol write batch_script` and preserved it as
`ob_replay_batch_command_20260916.sh`. Commit this amendment and command before
resubmitting normally without exclusive allocation/dependency. YGOB job
20917 and unrelated jobs were not modified. V2 simulation inference remains
active; its partial results have not been scored.

After pushing scheduling amendment `5dadc80`, resubmitted the recovered
batch command as **21088** with 32 CPUs, 64 GiB, two hours, shared-node
allocation and no dependency. `scontrol` confirms OverSubscribe=OK and null
dependency; `squeue` subsequently confirms it RUNNING. Existing pinned
launcher/core/input/output arguments are unchanged. Original 20919 remains
cancelled, not a failed or completed scientific run.

Added `prepare_orthobench_factorial.py` for the next step after replay success.
It gates on scheduler completion, equivalent stage partitions and provenance;
verifies the launcher's complete core source set against frozen 7f3a9e4;
checks cached gene/species classes against FASTAs; and prepares four separate
profile/candidate arms. Satellite expansion uses the production helper and
rebuilds/validates merge constraints independently for each profile arm.
Eight planned cells keep the prespecified tree/root/pair settings, with no
reference-scoring argument in inference commands. Profile-off still retains
the initial HMM search and must not be labeled HMM-free.

Four new tests cover factorial commands/own-arm constraints, cache species
ownership, independent expansion traces, preserved seed files, and rejecting
incomplete partitions. These and the existing phylogeny replay tests pass
(18 total). No candidate preparation or new factorial accuracy has run yet;
pin this builder and submit it behind successful replay 21088. Unconstrained
membership, matched sequence-search and QfO controls remain required.

## Native Profile Runtime Failure: Replay And Simulation Correction Required

Replay 21088 is terminal FAILED (1:0, 00:03:49), despite native inference exit
zero. Non-profile partitions match exactly; both profile partitions instead
equal their non-profile counterparts. All profile counters are zero. Direct
loader and isolated synthetic-cluster probes establish missing `pair_align.so`
in frozen checkout 7f3a9e4. The profile worker suppresses the resulting OSError.
See `PROFILE_RUNTIME_FAILURE_20260916.md` and the machine-readable probe.

This supersedes the proposed submission behind 21088 above: no factorial
preparation has run and that failed job cannot authorize it. Added a fail-fast
exact-checkout/interpreter profile smoke gate to the current replay launcher;
pinned historical launchers and active inference sources remain untouched.

Both simulation manifests use this checkout. All 140 fixed-length OrthoHMM
records and 96 available variable-length records at the audit have zero built
profiles. Marked fixed-length reports as defective-runtime diagnostics and
updated the claims checklist. Original machine-readable outcomes remain intact.
Do not score variable-length OrthoHMM as the intended frozen configuration.
Array 21010 remains active to complete reusable OrthoFinder results; no core
or runtime files were changed underneath it and no variable-length accuracy
was inspected. Repaired OrthoHMM runs and amended execution provenance remain
required. The independent constant-length OrthoFinder failure remains real.

YGOB 20917 was verified PENDING with elapsed zero; used pending-only scancel.
Accounting confirms CANCELLED by 1000, 00:00:00. Recovered batch script retained
as `ygob_cancelled_batch_20260916.sh`. YGOB scientific settings/reference remain
frozen and no accuracy has been inspected. Resubmit after separate-checkout
native build/provenance, runtime smoke, and corrected label-blind replay gates.

Previous full-suite handle 24582 was missing on recheck; no test result inferred
from that. New full-suite run found one test expectation too narrow for editable
import hooks (635 passed, one failed): wrong-checkout ValueError is also a valid
failure mode. Corrected that assertion while retaining exact-source rejection.

Final full suite: **636 passed in 22.98s**; scoped `git diff --check` clean.
Actual isolated probe rejects the incomplete frozen checkout with OSError and
passes the development checkout with the identical profile Python source hash,
recording `pair_align.so` hash and a 20-position synthetic profile. This positive
control does not replace the required separately built frozen runtime.

## Separate Frozen Native Runtime Built

Previous goal turn made substantive progress: diagnosed missing native profile
runtime, guarded replay, invalidated affected claims, and pushed 6b0f8e2.
Re-read the full goal and rechecked live jobs before continuing. Array 21010
remains active; original frozen source/runtime paths are unchanged.

Created separate detached checkout `benchmarks/work/publication_method_native_v2`
at the exact frozen 7f3a9e4 revision. Added a fail-closed CPU runtime builder
that refuses existing binaries/provenance, verifies clean frozen sources,
records compiler/source/binary hashes and commands, and exercises the native
profile stage. All three CPU kernels compiled successfully with GCC 13.3.0.
Synthetic profile construction passes with length 20; pair_align.so matches
the development control byte-for-byte. CUDA is explicitly absent.

Manifest: `publication_native_runtime_20260916.json`. New replay launcher
requires it and checks binary/source integrity both before and after inference.
Fourteen targeted runtime/probe tests pass, including missing/changed/extra
libraries, changed source, wrong revision and incomplete manifest rejection.
Amended the ablation protocol before new replay outcomes; preserve all v1
evidence and write v2 to a fresh output directory. No corrected simulation,
YGOB or factorial outcome has yet been produced by this build.

Full unit suite after runtime build/verification changes: **645 passed in
23.69s**. Scoped diff whitespace checks pass. Commit/push this runtime and
protocol before submitting the corrected label-blind replay from a pinned
launcher worktree. Simulation array 21010 remains active (tasks 63/64 running
at last live check); no variable-length accuracy has been inspected.

Pushed runtime/launcher/protocol milestone **36e45a2209fbc5e2fff644f66a7d08198dadfeff**.
Created detached launcher worktree `publication_native_launchers_v2` at that
revision and submitted corrected replay **21138**, confirmed RUNNING with
32 CPUs on bizon. Allocation remains shared, 64 GiB, two hours; no held-out
dependency. Recovered submitted batch script is
`ob_native_replay_v2_batch_20260916.sh`. It points to the pinned build manifest,
same historical audit and fresh `publication_ob_replay_check_v2` output.
No claim of equivalence until terminal verification; no factorial preparation
or corrected validation submitted behind an unverified outcome.

## Corrected Simulation Execution And Comparator Reuse

Previous goal turn made progress by building/pinning the native runtime and
launching replay 21138. Re-read the full objective and confirmed that job still
RUNNING before continuing. Original variable-length array 21010 is now fully
terminal: 67 tasks COMPLETED, three FAILED (60, 62, 69). This is scheduler
evidence only, not native admission or accuracy; no variable-length scores
have been inspected. Preserve all original outputs and executors.

Extended simulation preparation to require native build provenance and an
exact-checkout profile smoke. Optional comparator reuse verifies the original
manifest hash, identical generation/seed metadata, environments and scientific
arguments, then schedules only the two OrthoHMM modes. Reused comparator paths
remain unchanged. The executor checks only scheduled output paths for absence,
never launches reused tools, and records runtime verification before and after
each OrthoHMM process. Source-only OrthoHMM execution is now refused.

Extended the results assembler to reject source-only OrthoHMM admission and
require per-method native evidence for corrected runs. Mixed-provenance
assembly verifies both terminal arrays and their own pinned executors, reads
each run's original input/output/native evidence, and selects corrected
OrthoHMM plus original OrthoFinder rows. Scheduler IDs and method-manifest
hashes remain distinct; comparator failures are retained, not upgraded by reuse.

Prospective amendment `SIMULATION_RUNTIME_CORRECTION_PROTOCOL_20260916.md`
preserves both panels' scientific settings and statistical endpoints. The
corrected CPU build also enables compiled search kernels; that runtime change
is disclosed. No corrected simulation inference has been submitted. Replay
equivalence and newly pinned execution manifests remain prerequisites.

Targeted preparation/execution/assembly tests: 24 passed before mixed-row
support; final full suite after mixed-row support: **656 passed in 24.78s**.
Tests cover scientific-argument drift, changed generation/seed/comparator
settings, duplicate datasets, refusal to run without native evidence,
preservation of reused files, source-only admission rejection, mixed-row
provenance and mismatched truth rejection. No new accuracy was calculated.

Pushed **b66225dc5fc355702575aebe34b644916894f236** and pinned
`publication_native_simulation_v3` there. Materialized both corrected manifests
from that pinned builder, retaining 70 datasets each and scheduling only the
two OrthoHMM modes. Fixed manifest SHA256:
`523a603b3a3ba759d50ddd8cdcae540d148920d189de073bbf768c3d1b4c17ed`;
variable manifest SHA256:
`bf677728cf9312edd9755d9ac1ceefbce3abb8c3848d43413631dc00ca00eb2f`.
Pinned executor --check-only passes on actual baseline_20261001 and
baseline_20261101 inputs, including environment/runtime/history/input gates.

Replay **21138 completed 0:0 in 00:08:59**. All four partitions are byte-equal
and partition-equal to historical stages (no differences). Corrected runtime
built **15,901 profiles**, considered 5,014,697 candidates, retained 579,444
significant profile hits and 24,738 strict profile edges. Preserved lightweight
replay and verification snapshots as `ob_native_replay_20260916.json` and
`ob_native_replay_verification_20260916.json`. No reference accuracy was scored.
This satisfies the cached-stage equivalence gate, not the remaining publication
requirements or end-to-end search equivalence. Commit the two new simulation
manifests and replay evidence before submitting the selective rerun arrays.

After pushing manifest/equivalence milestone **62c18a2**, submitted corrected
variable-length OrthoHMM array **21142** and fixed-length array **21143**.
Both use the existing generic Slurm launcher in pinned executor
`publication_native_simulation_v3` (b66225dc5fc355702575aebe34b644916894f236),
the new manifest/hash above and their unchanged panel roots. Each array has
70 tasks, at most two concurrent, four CPUs and 16 GiB per task, 12-hour limit.
Only the two OrthoHMM modes execute; OrthoFinder outputs are reused separately.
Logs: `benchmarks/work/publication_native_simulation_logs_v3/{variable,fixed}_%A_%a.log`.
Live squeue confirms tasks 0/1 RUNNING in both arrays; remaining tasks pending
array limits. No corrected accuracy has been evaluated.

Next required work: native completion/scoring after both corrected panels are
terminal (reuse original arrays 21010 variable and 20957 fixed with their own
pinned executors), restore YGOB with runtime guards, and prepare/run the
OrthoBench factorial now that cached replay equivalence is established.
All remaining original publication requirements remain active.

## Factorial Preparation And YGOB Runtime Guard

Previous turn progressed through verified replay equivalence and selective
simulation submissions. Re-read the full objective and rechecked live arrays
21142/21143 before proceeding; neither corrected panel has been scored.

Submitted label-blind OrthoBench candidate preparation **21161** using pinned
builder b66225d, verified replay 21138, and fresh `publication_ob_factorial_v1`.
Job completed 0:0 in 00:01:03. Manifest is `prepared_not_reconciled`, with eight
planned cells. No-profile arm has 63,245 seed families and 54,745 expanded
candidates (8,500 merges); profile arm has 62,885 seed families and 54,445
candidates (8,440 merges). Each expansion produced and validated its own
membership trace. These are stage counts, not accuracy. Reconciliation cells,
unconstrained diagnostic, matched sequence control and QfO remain unfinished.

YGOB launcher now requires a pinned launcher revision and validated native
build manifest, verifies exact prepared file sets, and records synthetic
profile probes before/after each OrthoHMM stage. Fixed its ephemeral Slurm
spool-file provenance by preserving the exact submitted script in the output.
Records MAFFT/FastTree/DIAMOND/OrthoFinder entrypoint hashes and both Python
environments; confirms OrthoFinder 3.1.5 before inference. The label-blind YGOB
verifier checks the native manifest snapshot and all four exact-source/binary
probes in addition to existing completion/input/output checks.

Prospective YGOB amendment retains all scientific settings and references but
changes exclusive scheduling to shared 32-CPU/128-GiB/24-hour execution. Timings
are explicitly contended, not matched efficiency evidence. No YGOB outcomes
have been inspected and no existing result directory has been overwritten.

Validation: full unit suite **657 passed in 25.10s**; final YGOB verifier tests
**14 passed**; shell syntax and scoped diff checks clean. Pushed **db3b5aa**
and pinned `publication_ygob_native_v2` at full revision
db3b5aa9ff05eaa81da58b353d4bf6d426bab4bc. Its real-input --check-only verifies
the frozen inputs/source/native runtime and passes without creating inference
outputs. Frozen build manifest SHA256 remains
aebea83807356b02307473506fa30c2dbd2c511d7ba75a0655eb12180a474d74.

Submitted corrected YGOB **21192** with the amended shared-node allocation,
32 CPUs, 128 GiB and 24-hour limit. Live squeue confirms RUNNING; the first
compute-node `high_sensitivity_before.runtime.json` reports passed, exit zero,
the correct separate frozen checkout and profile length 20. Output remains
`benchmarks/results/ygob_validation_v1` (original job never started); log is
`benchmarks/work/ygob_native_validation_21192.log`. The job retains the same
prepared 16-species input/reference and scientific settings. It will execute
the two OrthoHMM modes and full OrthoFinder sequentially; all scoring remains
separate and gated. No YGOB or corrected simulation accuracy has been inspected.

Prepared OrthoBench candidate manifest snapshot is committed as
`orthobench_factorial_prepared_20260916.json`. Next: verify and execute the four
reconciliation cells from this manifest, then score all eight cells under the
prespecified protocol; retain required unconstrained and sequence controls and
the QfO factorial. Continue monitoring both corrected simulation arrays and
YGOB without changing active runtime/source files.

## OrthoBench Reconciliation Execution Wrapper

Previous turn progressed through candidate preparation and guarded YGOB
resubmission. Re-read the full objective and confirmed live simulation arrays
21142/21143 and YGOB 21192. None has been scored. Now implementing the four
R=1 cells from the prepared eight-cell factorial; the four R=0 cells use the
existing candidate partitions.

Added `run_orthobench_factorial_cell.py`: exact plan comparison rejects changed
CPU/settings, foreign constraints or reference-scoring arguments; verifies
preparation completion, replay equivalence, input/trace/source file hashes and
full core source sets; uses the pinned replay harness and frozen tool environment.
Reuses the tested execution recorder for refusal to overwrite outputs, per-cell
GNU-time/log/status evidence and complete output inventories. Dependencies are
rechecked after execution. Success remains pending native validation/scoring.
The execution scope is reconciliation only, with profile expansion inherited
from the verified upstream replay, not recomputed in a source-only checkout.

Seven targeted tests pass. Actual p0_c0_r1 --check-only passes against the full
OrthoBench inputs and frozen files, without inference or labels. Added a pinned
Slurm array launcher, four tasks with at most two concurrent, 32 CPUs/64 GiB
and 24 hours each. Prospective execution amendment records the shared-node
timing limitation before any reconciliation outcomes. Full tests pending.

Full suite **664 passed in 27.67s**; shell syntax and scoped diff checks pass.
Set the child working directory explicitly to the pinned replay worktree and
record it, so git provenance cannot inherit an unrelated submission directory.

Pushed execution milestone **9a8630197401d97e4fcc131f4939433bd3b89aac** and
pinned `publication_ob_factorial_execution_v1` there. Real expanded-cell
p1_c1_r1 preflight passes through the actual pinned Slurm launcher, in addition
to the previously checked unexpanded p0_c0_r1. Seven focused tests also pass
after explicit working-directory recording.

Submitted reconciliation array **21248**, four tasks with concurrency two.
Live squeue confirms tasks 0/1 RUNNING and 2/3 pending the array limit. Initial
cell status records confirm p0_c0_r1 running with the correct pinned replay
working directory `publication_native_simulation_v3`. Cell mapping is
0=p0_c0_r1, 1=p0_c1_r1, 2=p1_c0_r1, 3=p1_c1_r1. Outputs remain under
`publication_ob_factorial_v1/cells`, execution evidence under its `execution`
directory, and Slurm logs under `publication_ob_factorial_logs_v1`.

At submission, YGOB 21192 and corrected simulation arrays 21142/21143 were
still active. No ablation, YGOB or corrected simulation accuracy has been
inspected. Next: validate terminal outputs and frozen statistics, finish
matched sequence/unconstrained controls and QfO factorial, then continue the
remaining error, robustness, biological-application and publication work.

## Prespecified OrthoBench Factorial Statistics

Previous goal turn made progress by validating and launching reconciliation
array 21248. Re-read the full publication objective and confirmed active
simulation arrays, YGOB and factorial jobs before working on statistics.
No ablation, YGOB or corrected simulation accuracy has been inspected.

Added `bootstrap_orthobench_factorial.py`, consuming already-validated per-RefOG
sufficient statistics rather than reading inference outputs. It reuses the
audited weighted OrthoBench statistic and recomputes it for each shared RefOG
bootstrap draw. Defaults retain the prespecified 20,000 PCG64 multinomial
draws and seed 20260918. All eight cells are explicit; running/unknown cells
are rejected. Failed cells have reasons and no scores; their three affected
conditional effects are unavailable, with no zero imputation or reduction of
the 36-endpoint Bonferroni adjustment.

Implements all 12 conditional factor effects with nominal and adjusted F1/P/R
intervals and descriptive family wins/ties/losses. Adds the six descriptive
two-factor differences of differences, without interaction confidence intervals
or extra inferential claims, consistent with the original protocol. Renderer
includes all cells, all contrasts, failures and scope limitations. Profile-off
still retains HMM-based initial search; reconciliation changes output level.
Native completion/conversion/reference/provenance validation remains separate.

Initial ten tests pass, including independently recomputed scalar bootstrap
quantiles, weighted-vs-macro F1 distinction, zero effects for identical cells,
order invariance, unavailable contrasts without zero scores, and invalid family
or cell-state rejection. Added an explicit interaction-sign check and reporting
of zero generated draws when all cells fail. Full suite verification pending.

Full suite **675 passed in 29.54s**; final focused factorial tests **11 passed**
after the all-failed draw-count clarification. This statistical implementation
has not been applied to the running factorial. Corrected fixed-length array
21143 is no longer in squeue; next action is terminal-accounting-gated assembly
using pinned b66225d and original comparator execution evidence, not assuming
that absence from the queue alone proves scientific completion.

Pinned b66225d assembler completed the corrected fixed-length panel: all 280
explicit outcomes recorded in `simulation_fixed_native_results_20260916.json`.
Verified both corrected array 21143 and original comparator array 20957 with
their original executors and separate provenance. Generated the corresponding
table and `SIMULATION_FIXED_NATIVE_INTERPRETATION_20260916.md`.

High sensitivity: 70 admitted; satellite_v2: 64 admitted/six execution failures.
All original OrthoFinder outputs remain rejected by native gates; all 14
planned comparisons have zero complete pairs and no differences/intervals.
Every corrected OrthoHMM run built profiles (62-113), but none added profile
edges. All 280 status/score dictionaries match the defective-runtime snapshot.
Directly rechecked corrected failure logs: same missing connected single-copy
taxon coverage (n2 for seeds 3/6; n3/n4 for seed 7) in both divergent conditions.
No positive HMM-expansion or superiority claim follows from this panel.

Corrected variable array 21142 also left the queue. Started its separate
terminal-accounting-gated assembly using b66225d and original comparator array
21010/pinned executor f5f4e1b. Output is a new
`simulation_variable_native_results_20260916.json`; no partial scientific
summary has been inspected. Assembly session remains active at this update.

Variable native assembly completed: 280 explicit outcomes, separately verified
against original comparator array 21010/f5f4e1b and corrected array 21142/b66225d.
Generated `SIMULATION_VARIABLE_NATIVE_RESULTS_20260916.md` and interpretation
from `simulation_variable_native_results_20260916.json`. These are the first
inspected variable-panel accuracy results; no scientific settings changed.

High sensitivity admits 70/70, satellite_v2 67/70, full OrthoFinder and its
checkpoint each 65/70. Full OrthoFinder leads every paired condition mean.
Adjusted F1 intervals exclude zero below it for all seven high-sensitivity
contrasts and four satellite contrasts; satellite turnover/missing20/uneven_taxa
intervals include zero. Baseline satellite/full F1 is 99.46/99.94%; turnover
98.84/99.36%. Divergent satellite paired differences are -11.74 points (n=5)
and -12.28 for divergent_turnover (n=8), not differences of unmatched table
means. Recall is the main observed deficit, not a demonstrated causal mechanism.

Rechecked native failures: three satellite tree-coverage failures (seed 9 in
both divergent conditions, seed 10 in divergent_turnover); five OrthoFinder
nonfinite-graph failures in divergent seeds 1/2/7/8/9, with checkpoint exclusion.
Exact causes of those remaining comparator numerical failures need diagnosis.
All 140 OrthoHMM runs build 62-193 profiles but add no profile edges. No HMM
expansion advantage or overall-superiority claim follows. Both simulation
panels remain separate; YGOB and factorial inference/scoring and other original
publication requirements are still incomplete.

## Native Diagnosis Of Variable-Length Comparator Failures

Previous turn completed corrected simulation scoring and the prespecified
factorial statistics. Re-read the full objective and confirmed YGOB 21192
and factorial 21248 remain live before auditing the five variable-panel
OrthoFinder failures. No active inference sources or files were modified.

Added `diagnose_orthofinder_normalization.py` and ran it under the frozen
OrthoFinder interpreter on all ten divergent seeds, not only failures. Verified
native source/package provenance and all saved inference files before/after.
Recomputed 640 matrices using the native maximum-score BLAST reader with exact
self-hit exclusion and native length normalization. No scientific accuracy
was calculated or changed. All input proteomes have 82-92 distinct lengths.

Reproduced seven nonfinite within-species matrices in exactly the five failed
seeds; the five valid seeds produce none. Each affected matrix has two
non-self hits at one length product, giving rank-one two-parameter fitting.
Native intercept exponentiation overflows for the resulting extreme fits.
This is local degeneracy, not the old global equal-length assumption. Of 23
rank-deficient fitted subsets, seven are nonfinite, three have no stored
normalized scores, and thirteen are finite. No normalization call raises.

Result artifact `orthofinder_variable_normalization_diagnostic_20260916.json`
records all matrices, warnings, fitted parameters, source/input hashes and
graph correspondence. Interpretation in
`ORTHOFINDER_VARIABLE_NORMALIZATION_AUDIT_20260916.md` retains original
exclusions, acknowledges finite rank-deficient cases, and avoids silently
repairing or relabeling the competitor. Ten targeted tests and the full
**685-test suite pass** (32.84s); scoped whitespace checks pass. No upstream
issue was submitted. Native failure explanation is now supported, but general
failure frequency and remaining publication requirements are not established.

## Manuscript Simulation Integration And Factorial Admission Issue

Read the full objective and rechecked live scheduler state. The preceding
prompt-writing turn did not advance experimental state; this turn updates
the manuscript and diagnoses a newly terminal verification failure.

Added corrected simulation methods, runtime admission, paired-seed inference,
both panels' outcomes, normalization failures, and zero profile-expansion
edges to `PUBLICATION_MANUSCRIPT_DRAFT_20260916.md`. Removed the outdated
statement that multi-seed simulation remains wholly unfinished, without
removing the outstanding robustness, error tracing, controls or independent
validation requirements. Crosschecked reported percentages and contrasts
against the generated variable-panel results. Checked 26 single-line local
evidence links with Perl; multiline links also require review. Scoped
whitespace verification passes. This documentation change does not change
inference, scores, or frozen scientific settings.

Scheduler accounting: YGOB 21192 remains running; factorial 21248 tasks 1/2
run and task 3 waits. Task 0 is FAILED 1:0 after 26:37. Its preserved status
records `finished_pending_native_validation`, no failed methods, and native
process exit 0. The batch traceback is the post-execution environment check:
`Package inventory changed: orthofinder`, not an inference error.

Reproduced the exact inventory query under the frozen OrthoFinder interpreter:
from the repository root it matches the frozen inventory; from the b66225d
replay root it differs only by absence of `orthohmm: 0.5.0`. The factorial
runner changes its own cwd to the replay root before execution and does not
restore it before post-validation. This supports a cwd-dependent distribution
discovery defect, not an installed-package change. Preserve outputs and failed
scheduler accounting. Next: correct future verification context with tests,
perform independently recorded postflight recovery of all affected outputs,
then complete native validation/conversion before scoring. Do not silently
relabel the failed batch as successful or restart expensive inference.

## Factorial Runner Working-Directory Fix

Previous turn made progress by integrating simulation evidence and diagnosing
the first factorial postflight failure. Re-read the objective and confirmed
YGOB 21192 and remaining factorial tasks are live. Task 21248_1 has now also
terminated FAILED 1:0 (32:45) at the identical postflight inventory check;
its inference record reports exit 0, no failed methods, and 52,331 inventoried
artifacts. These are not yet admitted scientific results.

Changed only the development factorial runner to restore its original cwd in
a finally block after inference, before source/environment postflight checks.
Child execution retains the intended pinned replay cwd. Added tests asserting
the execution cwd, equal pre/post verification contexts, and cwd restoration
after successful inference, reported method failure, and raised exceptions.
All 10 targeted tests and the complete **688-test unit suite pass** (31.28s).
Scoped whitespace checks pass. The real prepared p0_c0_r1 check-only command
passes the frozen source, input, runtime and environment checks from the
original repository directory; no inference or scoring was performed.

Pinned running executors, original failed statuses, scientific settings and
outputs are unchanged. The corrected development runner does not retroactively
admit either failed batch. Required next work is separately recorded recovery
verification of finished artifacts and exact execution provenance, followed
by native completion, root-HOG conversion and all-cell scoring gates. A
successful check-only preflight is not a substitute for these output checks.

## Preserved Factorial Output Integrity Recovery

The preceding turn fixed and tested verification cwd restoration. Read the
full objective again and confirmed YGOB 21192 and factorial tasks 2/3 live.
Added `recover_factorial_postflight.py` for the specifically diagnosed array
21248 defect. It requires unique FAILED 1:0 scheduler records, exact task/raw
job identities, the known postflight traceback, successful method execution,
and unscored/unvalidated original statuses. Other failures and live tasks are
rejected. It does not change scheduler state or any inference artifact.

Rechecked the original 9a86301 executor revision and clean tracked benchmark
sources, recorded source/manifest provenance, exact cell command and inputs,
prepared source/FASTA/candidate checks, current native runtime/environment and
executable resolution from the original verification cwd. Verified every
recorded output hash and the complete output file set with `verify_process`;
original status and batch-log hashes remain unchanged during the audit.

Both completed cells pass integrity recovery: p0_c0_r1 has 56,021 artifacts;
p0_c1_r1 has 52,331. New separate records are
`orthobench_factorial_postflight_p0c0_20260916.json` and
`orthobench_factorial_postflight_p0c1_20260916.json`. Both retain FAILED
scheduler evidence and explicitly set native validation, accuracy evaluation,
and scoring admission false. No inference rerun or scientific score occurred.
These checks cannot retrospectively prove every transient execution state;
they establish preservation and present postflight integrity under the
documented cwd correction, not complete native-output semantics.

Twelve failure-gate tests and the full **700-test unit suite pass** (31.80s).
Actual-data audits cover both finished tasks; scoped whitespace checks pass.
Next: native completion/root-HOG conversion validation and all-eight-cell
scoring assembly, followed by prespecified factorial statistics. Remaining
tasks, YGOB independent validation, QfO ablations, matched search control,
robustness/scaling, biological application and publication packaging remain
unfinished. No full-goal completion claim is made.

## Native Root-HOG Partition Conversion

Previous turn completed preserved-artifact integrity audits. Read the objective
and verified YGOB 21192 and factorial tasks 2/3 remain live. Inspected native
replay metrics, reconciliation manifests, summary records and root-HOG writer
semantics. Implemented `validate_factorial_partition.py` using existing strict
native readers/membership checks and source-family split accounting.

The converter checks complete FASTA identifier coverage, duplicate membership,
native sequential HOG IDs, canonical candidate-family IDs, absence of
cross-source merges, agreement of candidate/root counts with native summary,
family completion counts, and replay/native summary equality. Singleton groups
are retained. Files are hashed before and after validation; converted groups
and a separate evidence record are generated without reading reference labels.
This checks partition semantics only, not every native reconciliation gate.

Both finished cells pass. p0_c0_r1: 63,245 candidates, 64,925 root HOGs,
676 split source families. p0_c1_r1: 54,745 candidates, 60,092 root HOGs,
2,026 split source families. Both preserve all 251,378 FASTA genes and have
zero cross-source merges. Evidence records are
`orthobench_partition_p0c0r1_20260916.json` and
`orthobench_partition_p0c1r1_20260916.json`; converted text remains under
`benchmarks/results/publication_ob_factorial_v1/validated_partitions/` rather
than source control. No inference outputs were changed and no accuracy was
computed. Full native source/tool/tree/membership checks and all-cell scoring
assembly are still required before admission.

The full collected unit suite passes **711 tests** (38.44s). An additional
valid-family-split case was added after collection; all **12 targeted converter
tests pass**, including that case. Scoped whitespace checks pass. These
structural counts do not establish an accuracy advantage or completion of the
ablation experiment or publication goal.

## Combined Native Group-Output Gates

Previous turn completed strict root-HOG conversion. Re-read the objective and
confirmed live YGOB 21192 and remaining factorial tasks 2/3. Added
`validate_factorial_native.py`, joining a fresh full artifact-integrity audit
with native replay command/source/input checks, exact scientific parameters,
summary agreement, candidate/constraint/output hashes, tool paths/versions,
inferred species-tree source/hash/taxon coverage and finite branch lengths,
membership accounting, and complete root-HOG partition validation.

Actual-data audits pass for p0_c0_r1 and p0_c1_r1. New records
`orthobench_native_p0c0r1_20260916.json` and
`orthobench_native_p0c1r1_20260916.json` retain the failed scheduler history
inside the successful postflight recovery evidence. Their scope is native
root-HOG group benchmark admission, not independent reconstruction of trees
or validation against pairwise truth. No reference labels or scores were read.

Expanded-cell native accounting matches all 8,500 supplied constraints:
5,928 supported and 2,572 detached. The unexpanded cell has no constraint
filtering. Both trees have complete native taxon coverage. Twelve targeted
metadata tests and the full **724-test suite pass** (36.42s); scoped whitespace
checks pass. Original executors, inference artifacts and scientific settings
are unchanged. Next: finish remaining cells' same gates and assemble all eight
prespecified OrthoBench cells, crosscheck scoring and run paired factorial
statistics. QfO ablations, independent validation and the other full publication
requirements remain outstanding; this is not publication completion.

## Eight-Cell Scoring Assembly Prepared

Previous turn validated native provenance for two finished cells. Read the
objective and rechecked YGOB 21192 and factorial tasks 2/3, all still live.
Implemented `assemble_orthobench_factorial.py`: it requires all four native
tasks terminal before reading outcomes, then freshly validates all four R=1
cells and the frozen R=0 candidate partitions before any scoring. A genuine
execution/integrity failure halts for diagnosis rather than silently dropping
a cell or imputing zero. No automatic inference restart is performed.

The assembly pins the previous reference snapshot SHA256 and official scorer
source, checks the complete 70-RefOG and 11-exclusion file sets, requires all
cells to partition the full FASTA universe, converts native groups explicitly,
and crosschecks weighted F1/P/R and exact-family counts against the native
OrthoBench CLI. The official CLI prints one decimal percentage place; its
rounding tolerance is 0.05000001 percentage points, not a precision claim.
It then uses the existing prespecified 20,000 paired-RefOG draws and 36-endpoint
Bonferroni factorial analysis. Inputs are rehashed after scoring and no result
directory may be overwritten. Genuine terminal inference failures currently
require a separately reviewed explicit failure path before analysis can proceed.

Fifteen new tests cover missing/duplicate/live task gates, rejection before
outcome reads, and official scoring agreement/rounding. Full **739-test suite
passes** (26.28s), with scoped whitespace checks passing. The actual assembly
command correctly exits at `Factorial still running; no partial scoring` and
creates no result directory. Independently checked the frozen reference and
official-source hashes without scoring predictions: 70 references and 11
exclusion files match. End-to-end scoring integration remains untested until
the final two cells are terminal and validated. No ablation accuracy has been
inspected and no completion claim is made.

## Third Cell Validated; Coverage And Resource Reporting

Previous turn implemented the all-cell scoring assembly. Read the full
objective and confirmed YGOB 21192 and factorial task 3 remain live. Task 2
terminated FAILED 1:0 at 25:30 with the same postflight verification defect.
Ran the unchanged combined native validator: p1_c0_r1 passes, with original
scheduler history preserved in `orthobench_native_p1c0r1_20260916.json`.
It retains all 251,378 genes across 64,616 root HOGs from 62,885 candidate
families, with 676 split sources and no cross-source merges. No accuracy read.

Extended the assembly report with per-cell coverage: all assigned genes,
singletons, nonsingleton membership, multispecies groups and their gene count.
The report explicitly distinguishes assignment including singletons from
orthology accuracy. Added native replay wall/user/system CPU measurements,
mean utilized CPU cores and sampled process-tree RSS, retaining byte units in
JSON and GiB only for display. Reject nonfinite/negative measurements and
incompatible RSS conventions. R=0 upstream partitions have null separate
resource measurements, never fabricated zero costs. All timings remain
incremental cached shared-node work, not matched end-to-end speed evidence.

Nine additional tests cover coverage categories, omitted genes, CPU arithmetic,
memory convention, invalid measurements and NA rendering. All 24 assembler
tests and the full **748-test suite pass** (33.02s). Scoped whitespace checks
pass. Fourth-cell validation and eight-cell scoring remain pending; no
inferential comparison is reported from the three available cells. YGOB
batch log now records the high-sensitivity metrics path, but the overall job
is still running and held-out outcomes remain uninspected.

## Corrected Simulation Evidence Figures

Previous turn validated the third cell and added coverage/resource reporting.
Read the full objective and confirmed factorial task 3 and YGOB 21192 remain
live. Used the completed corrected simulation results to generate publication
figures while these jobs run; no new inference or accuracy calculation.

Added `plot_simulation_evidence.py`, requiring an explicit source JSON hash.
Each panel shows native-admitted counts for all four methods and paired F1
effects for the two OrthoHMM modes, with included-seed counts and nominal and
Bonferroni-14 intervals. The fixed-length panel plots no zero-effect markers
where no comparator pairs are admitted. Captions disclose conditional success,
ten-seed uncertainty, nonpooled panels and parent-gated sequence checkpoints.

Generated PNG/PDF/SVG and manifests under
`figures_simulation_variable_native_v2_20260916` and
`figures_simulation_fixed_native_v2_20260916`, and linked both in the manuscript.
Viewed both initial renders, corrected a crowded repeated condition label,
then visually inspected both final PNGs: labels and intervals are visible
without overlapping neighboring panels. Initial drafts remain untracked and
are not the manuscript figures. Three plot tests pass, checking exact paired
interval coordinates, sample-size labels, no imputed effects with no pairs,
and rejection of unknown panels. These figures retain the negative simulation
findings and do not complete the remaining publication requirements.

## Development Profile Failures Now Surface

Previous turn completed and visually checked simulation figures. Read the
objective and verified final factorial task 21248_3 and YGOB 21192 still live.
Addressed the previously demonstrated silent native-library failure in the
development checkout only: `_build_profile_worker` now raises a contextual
RuntimeError, preserving the original exception cause and cluster ID, instead
of swallowing every unexpected exception as a missing profile.

The underlying builder's expected `None` cases remain unchanged, as do all
successful profile scores and scientific parameters. Four new regression
cases cover legitimate no-profile results, a missing-library OSError with
preserved cause, and unexpected failures in both serial and real spawn-worker
execution. Existing serial/parallel profile-equivalence tests pass. All
**13 profile-expansion tests** and the complete **755-test suite pass**
(34.12s). Scoped whitespace checks pass. Verified the frozen native benchmark
checkout's profile-expansion source remains identical to HEAD.

This release-oriented error-handling correction is not retroactively applied
to frozen benchmark results. Unexpected builder errors now stop inference;
expected no-profile returns can still occur and their causes/counts remain a
separate diagnostic concern. No active job, frozen checkout or scientific
output was changed. Remaining validation, ablation scoring and publication
requirements continue under the original objective.

## Scoring Assembly Integration Tests

Previous turn corrected development profile error handling. Read the objective
and confirmed final factorial task 21248_3 and YGOB 21192 still live. Added
two synthetic end-to-end assembly tests rather than reading partial benchmark
outcomes. They exercise all eight cell conversions, complete input membership,
all-four-native-validations-before-scoring order, crosscheck calls, 20,000
paired draws, twelve contrasts, JSON/Markdown serialization, coverage/resource
sections, and refusal to overwrite existing output.

Known fixture truth is a two-gene reference family: intact candidates score
100%, split native groups score 0%. The official-call boundary uses independently
specified expected values, not a second call to the same score function.
This is a harness integration test with mocked scheduler/native/reference
gates, not a claim that the real official executable or native inference was
tested by the fixture. A deliberate crosscheck mismatch leaves no published
JSON/Markdown result. Real-data official checks remain mandatory at assembly.

All **26 assembler tests** and the full **757-test unit suite pass** (30.63s);
scoped whitespace checks pass. The final factorial task remains RUNNING at
30:28 in the latest accounting poll. No benchmark results or scientific
settings were changed, and no partial factorial accuracy was inspected.

## Complete OrthoBench Factorial Scoring

Previous turn completed synthetic assembly integration tests. Read the full
objective and verified/waited on final task 21248_3 until terminal FAILED 1:0
at 32:55. Its traceback matches the known postflight cwd defect. Started the
all-cell assembly only after all tasks were terminal. All four fresh recovery
audits, native group-output validations and complete-gene partition gates pass.
Every cell matches official OrthoBench F1/P/R to printed precision and exact
RefOG counts. The 20,000 paired draws and 36-endpoint correction completed.
Assembly exited zero; no expensive inference was restarted.

First inspected complete factorial outcomes are now recorded in
`orthobench_factorial_results_20260916.json` and generated
`ORTHOBENCH_FACTORIAL_RESULTS_20260916.md`. Interpretation, manuscript and
claim checklist were updated. Reconciliation raises F1 in all four matched
settings with adjusted intervals above zero. Candidate expansion raises recall
and lowers precision; its adjusted F1 intervals include zero. Profile expansion
adds 0.568-0.704 observed F1 points, but all adjusted intervals include zero.
Full P1/C1/R1 reproduces historical F1 74.106074%; this is not a new OrthoFinder
superiority test. The initial HMM remains in profile-off cells.

All eight cells preserve 251,378 genes. Shared-node cached reconciliation takes
1,515-1,964 seconds with 1.467-1.577 GiB sampled summed process-tree RSS; no
matched end-to-end efficiency claim follows. Failed scheduler histories remain
explicitly retained despite no excluded scientific cells. YGOB remains live
and its outcomes uninspected. QfO factorial, matched search/unconstrained
controls, error tracing, robustness/scaling, biological application and the
remaining publication package are still required. Full objective remains active.

## QfO Replay Refinement Prerequisite

Previous turn completed the full OrthoBench factorial. Read the objective and
confirmed YGOB 21192 remains live. Began QfO ablation preparation by comparing
cached replay refinement with production `_refine_cluster_file`. Production
omits directed search-hit arrays at 50 or more dataset species, retaining
weighted graph edges for broad copy-only refinement. The replay supplied
directed arrays regardless of species count in both refinement stages.

Corrected the development replay to use the production threshold constant
and unique species count, preserving original arrays below the threshold and
empty lists above it. Both multipass and post-profile refinement use the
same selected arrays; initial HMM search, edge construction and profile
search inputs remain unchanged. Reports now record `refinement_directed_hits`.
The completed 12-species OrthoBench branch is unchanged, and its pinned
executor/output artifacts were not edited.

Six added tests cover 12/49/50/51/100 species with sparse repeated labels,
array identity below threshold, and many genes from only one species. Existing
replay CLI tests remain present. Full **763-test unit suite passes** (23.95s)
and scoped whitespace checks pass. This fixes a necessary replay mismatch;
it is not proof of QfO production equivalence. Next: inventory and validate
the actual QfO normalized-hit checkpoint, freeze a corrected replay launcher,
reproduce production outputs before preparing the QfO factorial. No QfO
ablation or held-out validation accuracy was inspected this turn.

## QfO Numeric Checkpoint Audit

Previous turn corrected the broad-panel replay branch. Read the objective and
confirmed YGOB 21192 remains live. Located the retained QfO checkpoint under
`qfo_benchmark/results/orthohmm_high_sensitivity_isolated/output/orthohmm_working_res/high_sensitivity_checkpoint`.
It is the production numeric checkpoint, not the pickle dictionary currently
accepted by the cached replay. Reuse requires a numeric input adapter rather
than reconstructing an 88-million-entry Python dictionary.

Added `audit_accuracy_checkpoint.py` with exact inventory, pinned manifest,
all-file SHA256 verification, mmap loading, bounded-chunk shape/dtype/index and
finite-score checks, unique gene IDs and descriptive self-hit accounting.
Ran against the historical manifest hash b90c787f...: all checks pass.
`qfo_numeric_checkpoint_audit_20260916.json` records 976,504 genes, 88,729,858
hits, 78 species and lexically sorted gene names. There are 976,195 self-hits,
zero nonpositive scores, and score range 0.0047840368-7.4837209302. Self-hits
are observed cached evidence, not automatically discarded or classified as
invalid. The eventual replay must preserve the same production treatment.

Seven targeted tests pass and scoped whitespace checks pass; the real-data
audit completed successfully without modifying checkpoint files or evaluating
orthology accuracy. This verifies numeric integrity, not FASTA/source matching,
hit completeness, duplicate ordered-pair absence or production replay
equivalence. Next: add and test a numeric-checkpoint replay adapter, freeze
the exact source/runtime/input manifest and run label-blind QfO equivalence.
No search rerun, QfO ablation scoring or held-out scoring was performed.

## Numeric Checkpoint Replay Input

Previous turn audited the retained QfO numeric checkpoint. Read the objective
and confirmed YGOB 21192 remains live. Added mutually exclusive
`--accuracy-checkpoint` and legacy `--hits-pickle` replay inputs, requiring
`--checkpoint-sha256` with numeric checkpoints. Numeric input runs the exact
inventory/hash/chunked-array audit, then retains the existing native gene
indices, species codes, score values and self-hits as read-only memory maps.
No large Python hit dictionary or numeric reindexing is introduced. Legacy
pickle indexing and its input provenance format remain available.

Five added regression cases cover unsorted native gene order, readonly mmap
arrays and self-hits, wrong manifest hash, missing/conflicting input modes,
and rejection of a missing hash before creating outputs. All 24 targeted
replay/checkpoint tests and the full **775-test suite pass** (25.45s).
Scoped whitespace checks pass. Actual QfO adapter integration also passes:
976,504 names, 88,729,858 hits, four readonly memmaps and 976,195 self-hits.

This completes numeric input adaptation, not graph/profile replay equivalence.
Before the QfO run, freeze a launcher using this adapter and the broad-panel
refinement correction while preserving the intended frozen core algorithm;
verify FASTAs, historical source/runtime and target partition provenance.
No inference search, graph replay, QfO scoring or held-out scoring was run.

## Historical QfO Source And FASTA Binding

Previous turn completed numeric replay input. Read the objective and confirmed
YGOB 21192 live. Added `audit_qfo_replay_inputs.py` and ran it against pinned
historical metrics fb6b8d7e..., numeric manifest b90c787f..., source revision
694a77fe56167754ca949751bca88aa7d11353dc and target partition 63ade2f3....
All 30 recorded source files match their original Git commit, despite the
historical harness's dirty-worktree flag. All 78 FASTA hashes match. Every
one of the 976,504 checkpoint genes occurs exactly once in those FASTAs,
and each proteome maps bijectively to one cached numeric species code.
The final historical target partition hash also matches.

`qfo_replay_inputs_audit_20260916.json` records these bindings and seven changed
recorded source files relative to the newer frozen checkout: accuracy,
argument processing, main pipeline, parser, refinement, profile expansion,
and benchmark_production. This compares the historical recorded file set;
it is not a claim that no new files or native-library differences exist.
The historical source audit covers the recorded Python files, not an
unrecorded historical binary/environment inventory.

Six mapping tests and the full **781-test unit suite pass** (24.95s), plus
actual-data audit and scoped whitespace checks. This provides concrete input
and target provenance for the replay freeze; seven changed recorded sources
mean equivalence cannot be inferred from checkpoint integrity alone. Next:
pin the corrected adapter/runtime while keeping intended core settings,
execute the label-blind cached QfO replay and compare its native partition.
No graph inference or benchmark accuracy was evaluated this turn.

### QfO isolated replay launcher freeze (2026-09-16)

The immediately preceding conversational turn supplied a goal prompt, not an
analysis milestone. This continuation reread the actual objective and confirmed
YGOB job 21192 live. Created the isolated branch
`publication/qfo-replay-native-v1`, commit
`49ab110358c0b4c73806a640de9068494a311f63`, based on the corrected numeric replay
adapter at 9effec3. Only in that dedicated checkout, restored the historical
profile-worker exception behavior so all core sources match the frozen 7f3a9e4
publication implementation. The main development branch retains its explicit
profile-error reporting fix. Do not merge this launcher branch into development.

Copied the three already-verified CPU native libraries into the isolated
checkout without rebuilding them. `verify_qfo_replay_launcher.py` checks pinned
Git revisions, tracked source cleanliness, complete source/native file sets,
byte equality, the original native runtime manifest, and an actual profile
construction probe in the launcher interpreter. The successful evidence is
`qfo_replay_launcher_20260916.json`. It explicitly retains the limitation of
historical exception-to-None behavior; the probe is not per-cluster validation.

The isolated replay/input audits pass 30 targeted tests. Six new verifier tests
exercise source changes, native changes, added files, missing libraries, and
the exact-match path. An initial test command used the wrong tests directory;
the corrected tests/unit command passed. No QfO inference or accuracy scoring
has been launched by this milestone. Next: use this immutable launcher in a
no-overwrite, pre/postflight-verified QfO batch, compare the final partition
against the audited historical target, and retain non-equivalence if observed.

Full development unit suite: **787 passed in 23.95s**. Isolated launcher branch
was pushed successfully. GitHub still reports 21 dependency vulnerabilities
(1 critical, 7 high, 11 moderate, 2 low); triage remains required and no frozen
runtime dependencies were modified.

### QfO replay launched; YGOB inference completed (2026-09-16)

Previous turn was progress: isolated launcher freeze and verification. This
turn reread the objective and implemented `run_qfo_publication_replay.py` with
no-overwrite outputs, pre/postflight frozen launcher checks, fresh historical
input audits, installed-package version consistency, and full-universe final
partition comparison. It does not read benchmark accuracy labels. Parameters
are CPU32, BLOSUM62, CPM0.1, Leiden seed4, one profile pass, minimum one species,
and no jackknife. GNU time measurements are explicitly incremental and shared
machine; max RSS is not a simultaneous process-tree memory measurement.

Pinned executor commit `1c23311f8b01e0ef8fbfe499614bf62054560f54` is checked out at
`benchmarks/work/publication_qfo_replay_executor_v1`. The core/replay launcher
remains isolated commit 49ab110. Submitted the committed
`qfo_native_replay_batch_20260916.sh` as Slurm **21288**, confirmed RUNNING,
32 CPUs, 128GiB, 24h limit, zero restarts. Output directory:
`benchmarks/results/publication_qfo_replay_check_v1`; batch log:
`benchmarks/work/qfo_native_replay_21288.log`. The fresh audit has verified all
976,504 genes and 78 FASTAs. Replay equivalence remains unproven until native
inference and postflight finish. No results have been scored or tuned.

YGOB job **21192** became COMPLETED 0:0 after **1:41:20**. Ran the existing
label-blind verifier successfully and recorded
`ygob_native_files_verified_20260916.json`. This verifies scheduler/native
completion, frozen input and source records, OrthoHMM output manifests,
four profile-runtime probes, native library identity, tool entrypoint hashes,
and OrthoFinder completion/input copies. `all_scoring_gates_verified` remains
false: exact commands/versions, native conversion, overlap/reference-resource
checks and independent reference reconstruction remain before held-out scores
can be inspected. Do not equate this file gate with completed validation.

Seven new replay command/partition/no-overwrite tests pass, and the full unit
suite passes **794 tests in 24.14s**. Next: monitor 21288 without restarting;
complete the remaining YGOB admission gates, then evaluate the frozen held-out
panel. QfO factorial preparation depends on the replay result. The full
publication objective, including matched-search controls, robustness/scaling,
error analysis, biological application and archival deliverables, remains open.

### Independent YGOB reference reconstruction (2026-09-16)

Previous turn was progress: QfO launch and YGOB file verification. Reread the
objective and confirmed QfO job 21288 still RUNNING (4:36 at latest check).
Read the frozen YGOB protocol and recorded launcher. The two native OrthoHMM
commands agree with the declared CPU32/eight-worker-thread, BLOSUM62,
E-value1e-4, Leiden/CPM0.1, high-sensitivity/default-refinement settings;
satellite additionally records the prescribed inference/rooting/pair rules.
OrthoFinder log explicitly records version3.1.5, 32 search and eight algorithm
threads, and full default MSA tree inference. These observations still need
integration into the final machine-checked command/conversion admission gate.

Added `verify_ygob_reference.py`, a second reconstruction that does not import
the original preparation code or pillar parser. It transcribes retained
columns from the acquired README, independently detects duplicated genes and
excludes their entire rows, applies the frozen FASTA filters, and compares
every prepared protein's sequence/species and every reference membership.
Pinned raw snapshot hashes are checked before reconstruction. Actual-data
verification passed: 83,404 proteins, 16 species, 83,391 reference genes,
10,250 groups, excluded rows113/9896 and 13 retained input-only genes.
Evidence: `ygob_reference_reconstruction_20260916.json`.

Six targeted tests cover column/species mapping, genus/OFF filtering,
terminal-stop normalization, whole-row duplicate exclusion with retained
inputs, invalid FASTAs and malformed column counts. This establishes an
independent implementation check, not independent biological ground truth.
No accuracy predictions or outcomes were read. Remaining YGOB gates include
machine-checked commands/versions, native conversions, overlap/resource
audits and the frozen score/uncertainty assembly. Do not retune on YGOB.

Full unit suite: **800 passed in 24.63s**; scoped whitespace checks passed.

### YGOB native commands and group conversion verified (2026-09-16)

Previous turn was progress: independent reference reconstruction. This turn
reread the objective and confirmed QfO 21288 live. Added the label-blind
`verify_ygob_native_outputs.py`: reruns native file and reference checks,
requires exact native and harness OrthoHMM commands, verifies the recorded
OrthoFinder3.1.5 package and executed GNU-time command, confirms full MSA tree
inference in its log, and validates the four frozen native group conversions.

The first actual-data check rejected the raw global-numeric MCL filename
because the existing converter requires species/sequence ID pairs. No score
or final evidence file was emitted. Corrected the selection to the documented
`clusters_OrthoFinder_I1.2.txt_id_pairs.txt`; added an explicit filename
regression test. This fixes the new admission script, not an existing result.

Actual validation now passes, recorded in `ygob_native_conversion_20260916.json`:
all four methods contain exactly **83,404 unique input genes**, no foreign IDs
and no missing genes. Group/singleton counts: OrthoHMM high-sensitivity
8,665/3,321; satellite_v2 10,079/3,598; full OrthoFinder 7,343/1,863;
sequence-only checkpoint 6,845/1,754. These counts are label-independent and
are not accuracy outcomes. Native files have not been modified.

Eight tests exercise commands, foreign/duplicate membership, explicit missing
coverage, unknown methods, ambiguous file discovery, and correct MCL variant.
The report deliberately keeps all_scoring_gates_verified false. Remaining:
finish the overlap/reference-resource admission audit, then assemble and
independently check the frozen scores and paired uncertainty. Dependency
inventories cover recorded packages/entrypoints, not every executable invoked
internally by third-party tools; shared-machine timing limitations persist.

Final full unit suite: **808 passed in 25.72s**; scoped whitespace checks pass.

### Frozen YGOB evaluation completed (2026-09-16)

Previous turn was progress: native command/conversion verification. Reread
the objective and confirmed QfO 21288 still RUNNING (15:20 at latest check).
Located the already-completed prospective reference-resource audit rather
than inventing a new independence criterion. Added
`verify_ygob_overlap_screen.py`: rechecks retained hit/mapping/input/source
hashes, DIAMOND command/version, terminal screen20918, reviewed biological
input sources, and recomputes the descriptive overlap summary. All checks
pass: 71,714/83,404 proteins and 6,952/10,250 pillars have qualifying hits.
The evidence permits bounded novel-taxon transfer, not family independence.
Historical intermediate database byte identity is not proven by this check.

Ran `assemble_ygob_validation.py` after fresh native and overlap admission.
The frozen four-method scorer and 20,000 paired pillar bootstrap completed;
independent explicit pair enumeration matches every TP/FP/FN count. Full
per-pillar results and admission evidence are preserved at
`benchmarks/results/ygob_frozen_evaluation_v1/`. The compact committed snapshot
`ygob_frozen_results_20260916.json` links the full 7.9MB report by path/hash.
Markdown and interpretation: `YGOB_FROZEN_RESULTS_20260916.md` and
`YGOB_FROZEN_INTERPRETATION_20260916.md`. Manuscript and claims updated.

First held-out outcomes: satellite_v2 F1 **92.233654%**, high sensitivity
**82.038608%**, full OrthoFinder **92.318524%**, sequence checkpoint
**85.572844%**. Primary satellite-minus-full difference **-0.084870pp**,
nominal interval [-0.487503,0.312676], Bonferroni-six interval
[-0.622528,0.445222]. No superiority or formal equivalence established.
Satellite precision is higher by6.401035pp and recall lower by6.914033pp;
both adjusted intervals exclude zero. High sensitivity is worse on all
three endpoints with adjusted intervals excluding zero. All methods have
100% reference-gene coverage. No parameters, exclusions or endpoints changed.
Further outcome-driven development requires new independent confirmation.

Nine focused enumeration/screen tests pass; full suite **810 passed in
25.70s**. This completes the bounded frozen validation experiment, not the
overall publication goal. QfO factorial, matched sequence-search controls,
robustness/scaling, error analysis, biological application and archive/release
requirements remain open; raw YGOB redistribution permissions unresolved.

### Frozen YGOB publication figure (2026-09-16)

Previous turn was progress: completed the bounded held-out evaluation. Reread
the objective and confirmed QfO replay21288 live. Added
`plot_ygob_validation.py`, consuming the hash-pinned frozen result snapshot
3927d81b1fd851eef3c8ce1343679558ef4f93bcfe2e643435e369587a44a4ef.
It requires the admitted four-method panel and frozen bootstrap specification.
The figure shows all12 observed F1/precision/recall scores and all6 prespecified
paired differences, displaying nominal and Bonferroni intervals separately.
Its caption discloses complete coverage, the diagnostic checkpoint, the
curated-group endpoint and non-family-disjoint scope. The primary F1 interval
including zero is not presented as equivalence or superiority.

Generated PNG/PDF/SVG plus source/result/output hash manifest in
`figures_ygob_frozen_20260916/`. Visually inspected the PNG: axes, labels,
numeric annotations, interval marks and footnotes render without clipping or
overlap. Linked the figure and caption in the manuscript draft. Four tests
check panel contents, admission rejection, multiplicity and nonfinite scores.
No analysis settings or scientific results were changed.

Full suite: **814 passed in 25.46s**. QfO21288 is RUNNING at19:03 and
scontrol reports its live batch/process IDs. Both multipass and refined
multipass partitions now exist; profile-stage completion and final partition
equivalence remain pending. Do not treat those intermediate files as final
successful inference.

### Unconstrained satellite diagnostic launched (2026-09-16)

Previous turn was progress: YGOB figure. Reread the objective, confirmed QfO
21288 live, and implemented the outstanding unconstrained satellite control
named in the prospective ablation protocol. The current factorial executor
has a narrowly scoped --unconstrained mode, restricted to p1_c1_r1, with
separate output/evidence paths and the same frozen candidate partition,
core/replay launcher, environment, CPU32 and reconciliation settings.
Detailed exploratory three-endpoint comparison specification is now recorded
in the protocol; it follows factorial outcome inspection and is not presented
as independent confirmation. No YGOB outcome selects a new configuration.

First executor freeze3b885bf/job**21289** failed terminal1:0 after5s, before
native inference, because removing the constraints path was insufficient:
the frozen replay requires explicit --unconstrained-membership when a sibling
merge trace exists. Preserved its status/logs and original executor worktree.
Corrected the invocation, tested the exact flag, and moved outputs to
`p1_c1_r1_unconstrained_v2`. Executor freeze**d52add7** at
`benchmarks/work/publication_ob_unconstrained_executor_v2` was submitted as
**21290**, CPU32/64GiB/24h shared node. Native/core settings are unchanged;
this is a diagnostic, not a ninth core factorial cell.

New evidence/output paths are under
`benchmarks/results/publication_ob_factorial_v1/execution/p1_c1_r1_unconstrained_v2`
and `cells/p1_c1_r1_unconstrained_v2`, respectively. Batch logs are
`benchmarks/work/ob_unconstrained_21290.log`. Do not reuse the failed21289
directories or modify the completed original cells. Native validation and
official scoring must follow terminal success; no diagnostic score exists yet.
Four added tests prove unchanged original cells, exact argument-only changes,
and rejection of other factorial cells. All14 focused executor tests pass.

Corrected job21290 confirmed RUNNING and has created native phylogeny/working
directories. QfO21288 confirmed RUNNING at24:34. Full suite on the corrected
invocation: **818 passed in31.30s**. The prior 818-test run also passed before
the first launch, illustrating why unit/preflight success alone did not prove
the frozen replay would accept an incomplete diagnostic CLI invocation.

### Sequence-search control prepared (2026-09-16)

Previous turn was progress: launched the unconstrained diagnostic. Reread
the objective and confirmed QfO21288/unconstrained21290 live. Inspected the
built-in engine/profile implementation: initial HMM scores use raw integer
scores divided once by sqrt(query_length*target_length), while phmmer mode
uses different normalization and cannot directly enter high-sensitivity's
in-memory multipass branch. A sequence-search adapter is therefore required.

Read official DIAMOND command documentation and froze
`SEQUENCE_SEARCH_CONTROL_PROTOCOL_20260916.md` before control outcomes.
This is exploratory after development/held-out outcome inspection, not new
independent confirmation or a new default. Use DIAMOND2.1.11 very-sensitive,
12 species-specific target databases, all251,378 frozen query proteins,
CPU32, E-value1e-4, BLOSUM62 gaps11/1, composition1/masking1, all targets and
one HSP per pair. Record raw scores and both lengths; normalize with the same
formula but do not claim equal score calibration or sensitivity. Retain
self hits/asymmetry. Compare profile-off graph control to p0_c0_r0; a post
search top100-per-query/target-species diagnostic does not reproduce the
HMM prefilter cap. Two variants times F1/P/R form six exploratory endpoints.

Added and ran `prepare_sequence_search_control.py` against pinned OrthoBench
FASTA records. Prepared combined queries and gene/species/length metadata at
`benchmarks/work/ob_sequence_search_control_v1`, plus12 exact search plans.
Committed compact manifest snapshot `ob_sequence_search_prepared_20260916.json`
SHA9edb1bb8232e114e12a6199412dd45a31d740c8a1dcba36fcf05d2c58b23095b.
Four tests cover explicit settings/raw-score output, exact sequence/ID/length
preservation, no overwrite, changed inputs and duplicate IDs. No control
search or accuracy scoring has been run. Next: guarded search execution,
strict numeric adapter, hit-coverage diagnostics and frozen downstream replay.
Latest scheduler check:21290 RUNNING4:52;21288 RUNNING28:58.

Full unit suite: **822 passed in38.55s**; scoped whitespace checks pass.

### Sequence-search execution launched (2026-09-16)

Previous turn was progress: frozen protocol and actual input preparation.
Reread the objective and confirmed both existing jobs live. Added
`run_sequence_search_control.py` to execute the hash-pinned12-target plan.
It verifies preparation/source/input/query/metadata/binary records and exact
commands; refuses non-pristine output/evidence directories; records each
database/search phase's argv, exit code, elapsed time, stdout and GNU-time
log; preserves failures; and repeats input checks after all searches.
Successful execution is labeled pending numeric validation, never scored.
The two phase tests cover both successful and failed subprocesses and refusal
to overwrite evidence. Actual prepared-plan check-only validation passed.

Executor freeze**861eabf** at
`benchmarks/work/publication_ob_sequence_executor_v1` was submitted via
`ob_sequence_search_batch_20260916.sh` as **21291**, CPU32/128GiB/24h shared
node. It executes12 species-target searches sequentially. Confirmed RUNNING
and native DIAMOND alignment progress in target00/search.log; the first
database phase has completed. No control accuracy or hit-coverage summary
exists yet. Root evidence is
`benchmarks/work/ob_sequence_search_control_v1/execution.json`; batch log is
`benchmarks/work/ob_sequence_search_21291.log`. Do not rerun into that root.

Latest scheduler check:21291 RUNNING0:16;21290 RUNNING8:15;
21288 RUNNING32:21. Native logs provide actual search-progress evidence,
not merely a submitted-job claim. Next implement/validate numeric conversion,
then frozen graph replay and descriptive sensitivity/coverage comparisons;
do not infer equal sensitivity from equal E-values or similar score scaling.

Full unit suite: **824 passed in32.11s**; scoped whitespace checks passed.

### Sequence-search numeric adapter implemented (2026-09-16)

Previous turn was progress: launched the frozen search. Reread the objective
and confirmed all three jobs live. Added `convert_sequence_search_control.py`:
full conversion requires terminal-success job21291 and all12 successful
recorded phase commands plus frozen input/output hashes. It checks every
seven-field hit for known IDs, exact sequence lengths, target-species ownership,
finite positive scores and the frozen E-value cutoff. Raw score is divided
once by sqrt(query_length*target_length); bit scores are validated but not
mistakenly used as raw scores. SQLite primary keys reject duplicate directed
pairs instead of silently aggregating them. A disk-backed ranking supplies
the prespecified top100-per-query/target-species diagnostic, raw-score order
with lexical target-ID ties, independently of the all-hit variant.

Writes chunked memory-mapped arrays and native numeric checkpoints, then
applies the existing numeric audit. Both preserve the full input gene universe
and directed pair semantics; self hits are not explicitly discarded.
Original TSVs and frozen inference code remain unchanged. Refuses existing
outputs and verifies hashes again after conversion. Eleven targeted tests
cover normalization, malformed/unknown/wrong-length/wrong-species/nonfinite
hits, significance rejection, duplicate pairs, deterministic per-species caps,
and a real written/audited checkpoint including self hits.

Read-only parser validation of completed target00 passed for **2,948,026 hit
rows**. This is not full-panel numeric admission and produces no accuracy
result. Full conversion has not been launched while the12-target search is
still active. Latest scheduler evidence:21291 RUNNING4:37,
21290 RUNNING12:36,21288 RUNNING36:42. Next run guarded conversion after
search completion and implement the frozen graph-only replay/coverage checks.

Full unit suite: **835 passed in33.23s**; scoped whitespace checks passed.

### Sequence control downstream jobs queued (2026-09-16)

Previous turn was progress: numeric adapter and real first-target validation.
Reread the objective; confirmed all three existing jobs live and four completed
sequence-search targets. Available disk12TB; no storage blocker. Created
conversion executor freeze**8b8b0d3** at
`benchmarks/work/publication_ob_sequence_conversion_v1`. Submitted
`ob_sequence_conversion_batch_20260916.sh` as **21292**, afterok21291,
CPU8/64GiB/24h, output `benchmarks/results/ob_sequence_numeric_v1`.
The converter independently rechecks terminal parent success, not just Slurm
dependency satisfaction. This is queued work, not completed conversion.

Implemented `run_sequence_graph_control.py`, using the existing frozen49ab110
replay without FASTA inputs to disable profile expansion. No core algorithm
change is introduced. It requires completed conversion21292, both admitted
variants, recorded source/execution hashes and exact native runtime; each
variant receives separate no-overwrite outputs. Replay retains CPU32,
BLOSUM62, Leiden CPM0.1/seed4 and production graph/refinement. Native postflight
requires exactly multipass/multipass_refined stages, no profile expansion,
251,378 genes/12 species, complete final partition and unchanged runtime/
checkpoint hashes. No reference labels or scoring are passed to inference.

Executor freeze**588d92f** at
`benchmarks/work/publication_ob_sequence_graph_v1` is queued as array
**21293_[0-1%1]**, afterok21292, CPU32/64GiB/24h per task, serial variants
all_hits/top100. Output `benchmarks/results/ob_sequence_graph_v1/{variant}`;
batch logs `benchmarks/work/ob_seq_graph_21293_{0,1}.log`. Six new tests cover
the profile-free command and rejection of unexpected stages/universes.
Full unit suite **841 passed in37.71s**; scoped whitespace checks passed.

Latest scheduler:21293 and21292 PENDING Dependency;21291 RUNNING9:05;
21290 RUNNING17:04;21288 RUNNING41:10. No new completed inference or accuracy
outcome is claimed. After completion, validate native provenance, compare
hit coverage/sensitivity and score both controls against frozen p0_c0_r0.

### Unconstrained native admission implemented (2026-09-16)

Previous turn was progress: dependent conversion/graph execution chain.
Reread the objective and confirmed the five job handles in their expected
running/dependency states. Added `validate_unconstrained_control.py`, requiring
terminal-success21290 and its exact executor d52add73035ebe3472e3e1cb95c38a9421b09dee,
parent/treatment provenance, source/input/environment/tool inventories and
all execution artifact hashes before native admission. It does not accept
the failed21289 attempt or any already-scored/mixed-method status.

Extracted the existing native cell validator for reuse, preserving the
factorial checks and adding an explicit diagnostic policy branch. The branch
requires acknowledged unconstrained mode, no membership filter/constraint
accounting, the same frozen candidate input and rules, complete gene partition,
inferred species-tree coverage/finite branches, and recorded tool/output
provenance. Original constrained cells still require their own constraints.

Twenty-one focused native/status tests pass, including nine new acceptance/
rejection cases. Revalidated the actual completed constrained p1_c1_r1 cell
with the refactored validator: successful, recorded at
`orthobench_p1c1_native_recheck_20260916.json`; no accuracy recomputed. Full
unit suite **850 passed in33.91s**. The still-running diagnostic has not been
admitted or scored. Its validation CLI is ready for terminal completion.

Latest jobs:21288 RUNNING44:58,21290 RUNNING20:52,21291 RUNNING12:53;
21292 and21293 remain dependency-pending. Eight sequence-search targets are
complete and target08 is running. Preserve all active outputs and wait for
authoritative terminal evidence before conversion/scoring claims.

### Unconstrained scoring prepared; sequence search completed (2026-09-16)

Previous response supplied a goal prompt only, so it did not advance analysis.
Reread the active objective and revalidated live scheduler handles before
continuing. Job21291 completed successfully (0:0), elapsed17:32; its dependent
numeric conversion21292 is now RUNNING. Graph controls21293 remain pending
that conversion. This is terminal execution evidence, not yet numeric or
accuracy admission. The shared-node elapsed time is not controlled timing.

Added `assemble_unconstrained_control.py` for the prespecified exploratory
p1_c1_r1_unconstrained_v2 minus p1_c1_r1 comparison. Both native validators
must pass before predictions or reference labels are loaded. The assembler
checks full gene coverage, frozen reference resources, official scorer
agreement, and unchanged scoring inputs. It reuses paired RefOG bootstrap
statistics with20,000 replicates, seed20260918, and three-endpoint Bonferroni
adjustment; output is explicitly not a ninth factorial cell, independent
confirmation, a default change, or a publication-readiness claim. Resource
records retain their incremental/shared-node scope.

Seven new tests cover validation-before-label access, overwrite refusal,
comparison inventory, complete report assembly, official disagreement,
input mutation, and incomplete coverage. Full unit suite **857 passed in28.91s**.
No diagnostic accuracy has been evaluated while inference remains live.
Latest scheduler check:21290 RUNNING29:20,21288 RUNNING53:26,
21292 RUNNING3:49,21293 dependency-pending. Next: native-admit and score
the diagnostic after terminal completion; audit conversion and graph controls;
complete the QfO replay equivalence check without interrupting its active job.

### Factorial figure and unconstrained diagnostic completed (2026-09-16)

Previous turn was progress: gated diagnostic scoring implementation and tests.
Reread objective and polled authoritative scheduler state. Added publication
factorial plotter, seven tests, and PNG/PDF/SVG with hash manifest. It renders
all eight cells and all36 contrast endpoints, with nominal and adjusted
intervals; validated identities, point differences, finite scores and nested
intervals. Visually inspected the PNG: no overlapping labels or omitted
endpoints. Integrated figure into manuscript and corrected stale independent-
validation status to completed YGOB novel-taxon transfer with family overlap.
Full unit suite **864 passed in23.94s**.

During this turn21290 completed0:0 in32:44. Ran the gated assembler end to end:
both native outputs/provenance passed, complete input coverage was checked,
frozen reference admission and official-score crosschecks passed, and
20,000 paired RefOG draws with three-endpoint Bonferroni adjustment completed.
Full source output: `benchmarks/results/ob_unconstrained_scoring_v1`.
Committed copies: `orthobench_unconstrained_results_20260916.json` and
`ORTHOBENCH_UNCONSTRAINED_RESULTS_20260916.md`.

Unconstrained F1/P/R72.614510/76.092405/69.440641%, versus constrained
74.106074/81.770454/67.755336%. Differences -1.492/-5.678/+1.685pp;
adjusted intervals[-5.873,0.931]/[-14.307,-0.794]/[0.013,3.804].
Family F1 wins/ties/losses5/55/10. No F1 superiority or established F1 loss;
precision-recall tradeoff does not justify removing constraints by default.
Exploratory development-exposed evidence, not a ninth factorial cell or a new
independent confirmation. Manuscript records results and limitations.

Latest scheduler:21288 RUNNING59:09,21292 RUNNING9:32,
21293 dependency-pending. No numeric-control or QfO equivalence outcome yet.
Next: finish conversion/graph-control admission and sensitivity diagnostics;
complete QfO replay then its component analysis; retain remaining robustness,
error tracing, controlled scaling, biological application and release scope.

### Label-free search coverage implementation (2026-09-16)

Previous turn was progress: completed diagnostic scoring and publication figure.
Reread objective and confirmed21288/21292 RUNNING,21293 dependency-pending.
Implemented `compare_search_hit_coverage.py` to compare frozen HMM hits with
all-hit and post-search top100 DIAMOND variants before accuracy interpretation.
Requires successful terminal conversion21292, frozen input plan and cache
hashes, admitted numeric checkpoints, exact gene order and matching species
ownership. No benchmark reference labels are loaded. Rechecks input hashes
after analysis and refuses output overwrite.

Reports directed and nonself intersections, set-specific recovered fractions,
Jaccard, self hits, reciprocal directed/unordered pairs, queries without hits
or cross-species hits, target coverage, species-direction hit/query counts,
and normalized-score/count quantiles. Rejects duplicate or invalid pairs.
Checks top100 is a subset of all hits. Explicitly distinguishes hit overlap
from ground-truth sensitivity and score quantiles from comparable significance
thresholds. No accuracy, matched-efficiency or general HMM benefit is inferred.

Twelve focused tests passed; full unit suite **876 passed in27.12s**.
Prepared one-CPU64GiB dependent batch for an isolated committed executor;
submission details will follow. Conversion still actively writes its SQLite
database; did not inspect that database or restart any active process.

Submitted coverage job**21294**, afterok21292, from detached frozen executor
`benchmarks/work/publication_ob_hit_coverage_v1` at**4880d5d**. One CPU,
64GiB, four-hour limit; output `benchmarks/results/ob_search_hit_coverage_v1.json`,
log `benchmarks/work/ob_hit_coverage_21294.log`. Batch recorded in
`ob_search_coverage_batch_20260916.sh`. This diagnostic is label-free and may
run alongside graph inference; its timing is not an efficiency comparison.
No completed coverage result claimed until terminal and artifact validation.

### QfO replay completed but is not equivalent (2026-09-16)

Final polling showed21288 FAILED1:0 after1:03:45. This is the wrapper's
deliberate equivalence failure, not a native inference crash: native exit0,
postflight runtime/input/package checks passed, accuracy remains unevaluated.
Observed390,657 groups versus historical390,817;9,414 observed-only and9,574
expected-only groups. Final partition hash3b3ee97a1caef1775e7a5f89ff385a72316298428cbf5ca04ae493ebaf5b4fcc
differs from historical63ade2f317c6343bd0d1e98af0fe0e299dd61530fb6094dec62d97dab30a4df3.
Full evidence preserved in original output directory; committed copy
`qfo_native_replay_nonequivalence_20260916.json`. Do not reuse historical
QfO accuracy as a current frozen-core baseline or force equivalence. Next
investigate stage partitions, algorithm/source differences and determinism;
preserve both outputs and do not blindly restart the completed replay.

### QfO drift localized before profiles; targeted experiment prepared (2026-09-16)

Previous turn was progress: coverage implementation/queue and preserved QfO
nonequivalence evidence. Reread objective, confirmed conversion21292 live and
graph/coverage jobs dependency-pending. Compared historical/frozen source paths
and recorded stage counts. Initial RBNH edge counts agree24,148,515, but first
singleton-assignment counts differ1,361,622 vs1,360,934, preceding profiles.
This rules out profile scoring alone as the cause of the earlier count drift;
does not prove initial graph equality or identify the cause.

Added `QFO_REPLAY_DRIFT_DIAGNOSIS_20260916.md` with all recorded stage counts,
code-path findings and explicit limitations. Prepared targeted diagnostic
loading old/frozen RBNH implementations from exact Git blobs, checking graph
array byte equality and then repeating initial clustering twice on the same
preserved graph. Uses fixedCPM0.1/seed4, identical verified graph helpers,
full input universe, preserved partitions and singleton-array fingerprints.
No search, profile expansion, accuracy scoring or parameter tuning is rerun.
Historical intermediate partitions/native package inventory remain missing;
current repeatability alone cannot certify historical runtime equivalence.

Three new tests pass (one-bit weight/name mismatch, dtype/layout evidence,
actual old/frozen builder behavior on ties/self hits); full unit suite
**879 passed in26.13s**. One-CPU64GiB four-hour diagnostic batch prepared for
a committed frozen executor. Latest conversion check RUNNING21:40; no control
scores or hit-overlap results yet. Submission details follow separately.

Submitted initial-graph diagnostic**21295**, executor**1ff3154** at
`benchmarks/work/publication_qfo_graph_diagnostic_v1`, one CPU64GiB/four hours.
Output `benchmarks/results/qfo_initial_graph_diagnostic_v1`; log
`benchmarks/work/qfo_graph_diag_21295.log`. This is a targeted new experiment,
not a restart or replacement of failed-equivalence21288. Original outputs
remain unchanged. No graph-equivalence or repeatability outcome claimed yet.

### Sequence conversion completed; graph-control admission prepared (2026-09-16)

Previous turn was progress: targeted QfO diagnostic implementation and launch.
Reread objective and polled all active handles. Conversion21292 completed0:0
after25:12. The producer's final manifest records numeric checkpoints verified:
all_hits100,099,147 rows;top10048,991,663 rows. Both retain251,378 genes across
12 species and250,980 self hits; all scores positive and finite. Preserved an
exact committed manifest copy at `ob_sequence_numeric_conversion_20260916.json`.
These are search/conversion counts, not biological sensitivity or accuracy.
Both downstream chains started:21293_0 graph inference and21294 hit coverage;
the second serial graph variant remains dependency-pending.

Added `validate_sequence_graph_control.py`, requiring both terminal-success
array tasks with correct raw job identity, pinned executor588d92f, unchanged
source/runtime/commands/environment, admitted checkpoint contents and complete
gene partitions. Checks profile-off stage inventory, every stage hash/count,
final native-output hashes, and preflight/postflight provenance. Revalidates
artifacts after reading. Native admission does not evaluate benchmark labels.
Sixteen focused scheduler/status/command/profile/environment rejection tests
pass; full unit suite **895 passed in27.97s**.

QfO diagnostic21295 remains live with preserved RBNH arrays and its first
clustering repeat in progress. Do not infer final repeatability from these
intermediate files. Next: terminal admission of coverage and graph controls,
paired control scoring with six-endpoint correction, and completed QfO
diagnostic interpretation before choosing any additional replay experiment.

### Both sequence graph controls admitted; paired scoring prepared (2026-09-16)

Previous turn was progress: native validator and completed numeric conversion.
Reread objective and confirmed both21293 array tasks now COMPLETED0:0:
all_hits1:55,top1001:22, shared-node incremental replay times only. Ran
`validate_sequence_graph_control.py` on actual outputs; admission passed for
both variants, including checkpoint audits, source/runtime/command identity,
stage hashes/counts and full gene partition coverage. Evidence committed at
`ob_sequence_graph_native_validation_20260916.json`. No accuracy evaluated.

Added `assemble_sequence_search_control.py`: fresh native admission plus
successful terminal coverage21294 and its matching input/source/checkpoint
provenance are required before reference labels. Frozen HMM p0_c0_r0 candidate
partition and FASTAs are verified before comparison. Official OrthoBench
crosschecks and20,000 paired RefOG draws, seed20260918, six F1/P/R endpoints
with Bonferroni adjustment; no partial variant set or output overwrite.
Reports label-free hit diagnostics, gene coverage and graph resources with
process-only RSS/shared-node/incremental scope. Baseline graph cost remains
unmeasured, not zero. Does not assert equal E-values imply matched sensitivity.

Full unit suite901 passed in27.61s before adding resource extraction; all
eight focused assembler tests passed afterward, including two new resource
scope/invalid-value tests. Coverage21294 remains live5:29; QfO graph
diagnostic21295 live7:46. Prepared dependent scoring batch to recheck gates
when coverage completes; frozen executor/submission details follow.

Submitted scoring job**21297**, afterok21294, one CPU32GiB/four hours,
from frozen executor**b14b057** at
`benchmarks/work/publication_ob_search_scoring_v1`. Output
`benchmarks/results/ob_sequence_search_scoring_v1`; log
`benchmarks/work/ob_search_score_21297.log`. The batch starts from the original
repository verification directory and revalidates both graph tasks itself.
No score or coverage outcome claimed from a pending/live job.

### Sequence controls scored; tree perturbation panel prepared (2026-09-16)

Previous turn was progress: paired scorer/native admission and queued job.
Reread objective; coverage21294 completed0:0 in8:31 and scorer21297 completed
0:0 in0:30. Scorer repeated native/coverage/input/reference gates and official
crosschecks before20,000 paired RefOG draws/six-endpoint correction. Committed
exact result copies `ob_sequence_search_results_20260916.json` (includes full
coverage report) and `OB_SEQUENCE_SEARCH_RESULTS_20260916.md`.

HMM baseline F1/P/R69.763388/78.868592/62.542944%;DIAMOND all_hits65.762183/
55.036231/81.680884%;top10066.484803/55.897861/82.019046%. Control-minus-HMM
F1 -4.001/-3.279pp;adjusted CIs[-14.834,6.688]/[-14.424,7.844],both include0.
Precision decreases and recall increases with adjusted intervals excluding0.
Both controls family F1 wins/ties/losses39/10/21. No HMM-specific F1 advantage
established. Hit coverage differs substantially:18,235,373 HMM nonself hits
versus99,848,167/48,740,683; all-hit intersection covers65.12% of HMM hits.
HMM cache has no self hits, DIAMOND250,980, so raw no-hit query counts are not
commensurate. Cross-species/no-nonself coverage is retained for interpretation.
Manuscript and claims updated; no independent confirmation or default change.

Started separate robustness requirement: added prespecified rooted-NNI
protocol and deterministic topology generator using Bio.Phylo. Validated
actual p1_c1_r1 native tree, then generated unchanged supplied control, three
one-NNI and three two-NNI trees, selected by canonical topology hashes only.
Preserved taxa, branch-value multiset, explicit rooting and rooted clade
distances0/2/4. Roundtrip checks passed; source unchanged. Seven tests cover
neighbor enumeration, deterministic panel, topology/length preservation and
invalid inputs. Full unit suite **910 passed in27.25s**.
Prepared output `benchmarks/results/ob_species_tree_robustness_v1`; committed
manifest `ob_species_tree_robustness_prepared_20260916.json`. No reconciliation
or robustness accuracy has run; first require supplied-control equivalence
and validate checkpoint reuse. QfO graph diagnostic21295 still RUNNING13:12.

### Supplied-tree baseline replay prepared (2026-09-16)

Previous turn was progress: completed search controls and generated tree panel.
Reread objective; QfO diagnostic21295 remains live. Audited frozen pipeline
checkpoint handling: validated raw gene trees may be reused, but rooting and
reconciliation always rerun. Species-tree, checkpoint and output writes use
atomic replacement, protecting the original run when the replay helper seeds
hard links. Strengthened existing tree-change test to require unchanged raw
gene-tree hash and changed species-tree hash in the new reconciliation
checkpoint; passed alongside seven new command-treatment guard tests.
Full unit suite **917 passed in26.15s**.

Added `run_species_tree_control.py` for only the unchanged supplied-tree arm,
not the six perturbed trees. It verifies frozen tree/FASTA/candidate/environment
manifests and original native outputs before execution; retains constraints,
CPU32 and reconciliation rules; changes only supplied-tree mode/input, copied
checkpoint source and output destinations. Revalidates the original native
run after execution to detect mutation. Native equivalence and any scoring
remain separate gates. Output reserved at
`benchmarks/results/ob_supplied_tree_control_v1`; batch prepared for a pinned
executor with32 CPUs64GiB/four hours. Submission details follow. This cached
control is not end-to-end timing or completed tree-robustness evidence.

Submitted supplied-tree control**21298** from frozen executor**d9ea049** at
`benchmarks/work/publication_ob_tree_control_v1`,32 CPUs64GiB/four hours.
Log `benchmarks/work/ob_tree_control_21298.log`; evidence and new native output
under `benchmarks/results/ob_supplied_tree_control_v1`. No perturbation cells
launched: first require this arm's native validation and partition equivalence
to inferred p1_c1_r1. QfO graph diagnostic21295 remains a separate live job.

### Supplied-tree control equivalence established (2026-09-16)

Previous turn was progress: checkpoint-reuse test strengthening and supplied
control launch. Reread objective and revalidated live handles. Job21298
completed0:0 in2:35. Added `validate_species_tree_control.py`, requiring
terminal scheduler success, exact frozen executor/tree/candidate/environment
provenance, successful unscored execution/postflight and unchanged original
native checkpoint source. It reuses shared native checks, extended only for
explicitly declared supplied tree/checkpoint parameters and argv destinations.
Ordinary factorial cells retain inferred-tree requirements.

Actual full validation passed, including original p1_c1_r1 validation after
the shared-validator refactor. Supplied and inferred root-HOG partitions are
exactly equivalent:59,770 groups each, zero expected-only/observed-only groups.
Every reconciled family reused its raw gene-tree checkpoint:8,681/8,681;
rooting/reconciliation was recomputed. No benchmark accuracy recomputed.
Evidence: `ob_supplied_tree_control_validation_20260916.json`.

Eight new tests cover supplied metadata and execution-status admission;
full unit suite **925 passed in26.26s**. This establishes supplied-mode
baseline compatibility, not robustness to tree errors. Next freeze and launch
the six already-selected perturbation commands against this validated control,
then native-validate and score all18 prespecified F1/P/R contrasts. QfO graph
diagnostic21295 still RUNNING23:06 in its second clustering repeat; no final
repeatability result claimed from intermediate artifacts.

### QfO graph diagnostic completed; six tree perturbations ready (2026-09-16)

Previous turn was progress: supplied-control native admission and equivalence.
Reread objective. QfO21295 completed0:0 in24:18. Both old/frozen builders in
the current environment yield byte-identical24,148,515-edge graphs. Two
initial Leiden repeats yield identical349,898-group partitions and1,361,622
singleton edges, matching the historical count, not replay21288's1,360,934.
Independently verified recorded source/checkpoint/graph-file hashes, array
fingerprints and complete partition equality. Preserved exact report at
`qfo_initial_graph_diagnostic_20260916.json` and updated drift diagnosis.
This excludes the threshold-factor source change for this graph but does not
resolve replay drift. No archived historical initial graph/partition exists.
Next capture the exact replay entry point's initial graph/partition before
profile execution; do not claim general Leiden nondeterminism or historical
reproduction from the present evidence.

Extended the tree runner with six fixed perturbation indices. Each task
revalidates supplied-control equivalence, tree manifest, original native
source, inputs/environment, then retains original candidate constraints and
rules with only its tree/output changed. Reuses original raw gene trees,
recomputes reconciliation, and verifies source preservation afterward. Separate
outputs at `benchmarks/results/ob_species_tree_perturbations_v1/{nni1_0..nni2_2}`;
prepared array0-5%2,32 CPUs64GiB per task, four-hour limit. All18 planned
accuracy endpoints remain prespecified; none evaluated yet.
Eleven new index/identity/distance/command tests passed; full unit suite
**936 passed in26.24s**. Frozen executor/submission details follow.

Submitted tree-perturbation array**21299_0..5%2**, frozen executor**cfe0a09**
at `benchmarks/work/publication_ob_tree_perturbations_v1`. Logs
`benchmarks/work/ob_tree_perturb_21299_{0..5}.log`;32 CPUs64GiB/four hours per
task. Each task performs its own supplied-control and source-integrity gates.
No perturbation score claimed until terminal/native admission of the panel.

### Exact QfO replay initial-stage capture (2026-09-16)

The preceding prompt-only response made no analysis progress. Reread the full
objective and revalidated current state: tree array21299 tasks0-3 completed0:0,
tasks4-5 still live. No accuracy outcomes inspected.

Added a fresh-process observer of the frozen49ab replay entry point. It uses
the original cached replay arguments and environment overrides, records the
initial RBNH arrays and partition, calls the unchanged singleton builder, and
stops before the second clustering call. It compares the captured graph,
partition and singleton arrays with diagnostic21295. Before/after runtime,
input and installed-package checks guard interpretation. No profile searches,
refinement or accuracy evaluation run. Instrumentation itself is documented
as a limitation; this experiment is not a complete historical replay.

Five focused recorder tests passed; full unit suite941 passed in26.56s before
the final input/package guards, followed by five focused tests passing again.
Capture executor and Slurm submission will be recorded after freezing.

Capture frozen at**f6da2bc**, worktree
`benchmarks/work/publication_qfo_replay_capture_v1`, submitted as**21305**
(32 CPUs64GiB/four-hour limit). Confirmed live; historical input audit passed
976,504 genes/78 FASTAs. No capture comparison is available yet.

Tree array21299 now all COMPLETED0:0: tasks0-5 elapsed3:04,3:14,3:15,3:22,
3:50,3:33. Added an all-six native-admission validator checking exact frozen
executor, manifests, task/raw-job identities, unchanged supplied-control
equivalence, native output semantics and actual supplied/output topology.
Eighteen focused tests passed, including failure/live/identity/topology guards.
Actual panel admission is in progress; accuracy remains uninspected.

Panel native admission completed successfully: all six preserve251,378 genes,
54,445 candidate families and zero cross-source merges. Actual supplied and
native-output trees match their prespecified rooted clade distances2/4.
Root-HOG counts in fixed order are59,867;59,827;58,940;59,908;60,081;59,938.
These are output-integrity observations, not accuracy results. Exact admission
report: `ob_species_tree_perturbations_native_validation_20260916.json`.
Full unit suite **959 passed in25.11s**. Next implement and validate the frozen
18-endpoint official-score/paired-RefOG bootstrap comparison against supplied
control. QfO21305 confirmed RUNNING5:17; no diagnostic outcome yet.

### Tree robustness scoring prepared (2026-09-16)

Previous turn made progress: six native admissions, tested validator, and
frozen QfO capture submission. Reread full objective. QfO21305 remains live;
no restart or diagnostic result inferred from incomplete files.

Added a separately frozen scoring assembler for all six perturbations against
the unchanged supplied-tree control. Fresh all-panel native admission precedes
reference access; full input coverage and official-score agreement are
required. Uses the existing weighted70-RefOG statistic, paired20,000-replicate
bootstrap seed20260918, and all18 F1/P/R endpoints for Bonferroni adjustment.
Preserves family wins/ties/losses and incremental process-tree resource data;
does not select a best tree or alter defaults. Eleven focused tests passed.
Full unit suite and frozen submission details follow.

Full unit suite **970 passed in26.52s**. Prepared single-CPU16GiB/one-hour
scoring job with no output overwrite and frozen executor worktree.

Scoring executor frozen**1891846** at
`benchmarks/work/publication_ob_tree_scoring_v1`; submitted**21306**,
confirmed RUNNING1:23. Output `benchmarks/results/ob_species_tree_robustness_scoring_v1`.

QfO capture21305 completed0:0 in9:42. Independently reverified array hashes,
gene-name ordering, source/input records and complete partition comparison.
Identical24,148,515-edge RBNH arrays yield a different first partition:
349,950 groups versus diagnostic349,898 (4,427/4,479 unmatched groups).
Captured singleton edges1,369,532 differ from historical/diagnostic1,361,622
and earlier replay1,360,934. This localizes observed divergence to initial
clustering execution, before profiles; specific cause remains unproven.
Preserved capture/comparison JSONs and updated drift diagnosis. Next bounded
saved-graph worker repeats should record native binary/runtime identity;
no best-repeat selection or full-profile restart justified yet.

Tree scoring21306 completed0:0 in4:10. All seven native/scoring crosschecks
passed. Preserved exact result snapshot SHA256
8a00c80bb638c46ec37050f19caea3ed35cbf7488ef03578abead06d4d8b3412 at
`ob_species_tree_robustness_results_20260916.json` with generated Markdown.
Rechecked recorded assembler/reference/prediction hashes and all official
score comparisons. F1 spans73.744102-74.270946% versus supplied control74.106074%;
all18 adjusted intervals include zero. Family F1 ties62-69/70 per variant.
No equivalence/arbitrary-tree robustness claim and no new tree/default selected.

Added the all-variant score/18-effect figure with input hash verification,
seven malformed-evidence/rendering tests, PNG/PDF/SVG and provenance manifest.
Visually checked the rendered PNG: complete labels, no overlaps, all effects
visible. Full unit suite **977 passed in28.15s**. Integrated findings, figure,
QfO drift limitation and remaining requirements into manuscript/claim checklist.
Next QfO saved-graph worker identity/repeatability experiments, followed by
QfO component evaluation once reproducibility is understood; parameter/error/
application/scaling/release work remains required. Goal not complete.

### Preserved-graph worker repeatability diagnostic (2026-09-16)

Previous turn made progress: complete tree scoring/figure/manuscript and
localization of QfO initial-clustering divergence. Reread full objective;
jobs21305/21306 are confirmed terminal, not restarted.

Prepared three sequential fresh workers on the exact preserved QfO graph
from21305, without graph rebuilding, profiles or labels. Hash-pinned capture
records, source arrays, gene order and reference partitions are checked.
Each worker records imported graph/scientific module files, loaded shared
library hashes from/proc/self/maps, Python/package identity, relevant
environment overrides, CPU affinity and actual metadata. After native exit,
the parent rechecks files and complete partition coverage, comparing each
partition with both previous experiments and with the first new repeat.
All repeats remain in the record; failures stop without automatic retry.
Instrumentation and lack of historical binary identity remain explicit limits.

Eleven focused tests passed, including an actual fresh-process worker with
native libraries, expected exit behavior and an isolated gene. Full unit
suite and frozen submission details follow. Planned one CPU64GiB/one hour;
no end-to-end performance claim from this shared-node diagnostic.

Full unit suite **988 passed in25.33s** before freezing the repeat executor.

Frozen executor**a69d505** at
`benchmarks/work/publication_qfo_saved_graph_repeats_v1`; submitted**21307**.
Output `benchmarks/results/qfo_saved_graph_repeats_v1`, log
`benchmarks/work/qfo_saved_graph_repeats_21307.log`. No repeat outcome claimed
before terminal execution and source/partition validation.

### Prediction-independent error-feature preparation (2026-09-16)

Previous turn made progress: tested/frozen QfO worker-repeat diagnostic and
submitted21307. Reread full objective;21307 confirmed live, no restart.
While it runs, audited older `analyze_phylogeny_changes.py`: useful fixed
size/copy/identity concepts, but its pair statistics are not the publication
weighted OrthoBench statistic and cannot be substituted without validation.

Added sequence-only feature preparation against the hash-pinned12 OrthoBench
FASTA inputs, requiring all251,378 unique proteins. Explicitly records canonical
and noncanonical content, gaps/stops, length, global normalized entropy and
missingness. Shortness/composition descriptors are not biological fragment or
domain calls. No prediction/reference files or accuracy outcomes are read.
Twelve focused tests passed; frozen batch and full-suite evidence follow.

`ORTHOBENCH_ERROR_ANALYSIS_PROTOCOL_20260916.md` specifies the later exploratory
14-stratum/two-comparator/three-metric panel (84 adjusted endpoints), small-bin
and missing-data rules, reference-alignment validation, and all70-family tracing.
This is prospective for these joins but post-development, not independent
confirmation. Real duplication history, annotated fragments/domain architecture,
QfO extension, and biological application remain distinct unmet requirements.

Full unit suite **1000 passed in27.41s**, followed by12 focused tests passing
after adding scheduler identity and explicit zero-valued missingness counters.

Feature executor frozen**8bccb33** at
`benchmarks/work/publication_ob_sequence_features_v1`, job**21308** completed0:0
in23s. Verified source/input/table hashes,251,378 unique rows and exact
per-proteome counts. Preserved summary `ob_sequence_features_prepared_20260916.json`;
32.5MB TSV remains outside Git, SHA256
d315da458240fe2020b128b606adfcfdf15d4de34ee937a6b2aa7735f7452f97.
Input-only counts:13,822 short;2,345 composition-concentrated;460 composition
unevaluable;3,170 with noncanonical letters;177 with stops;no empty/gap/other
symbol sequences. No accuracy inference from these descriptors.

Located legacy70-RefOG alignments at
`benchmarks/results/orthobench_refog_alignments_20260902`; old generator uses
MAFFT --auto --thread1 but reuses any existing nonempty file without validating
IDs/sequences. Its identity statistic counts identical ambiguous symbols and
skips incomparable pairs. Do not reuse those summary values; current protocol
requires canonical-only comparisons and explicit missingness plus sequence/
provenance admission or a pinned rebuild. Family strata and outcome scoring
remain pending. QfO21307 remains running; no repeatability outcome claimed.

### Explicit reference-alignment preparation (2026-09-16)

Previous turn made progress: frozen error-analysis protocol and verified
sequence-feature inventory. Reread full objective;QfO21307 confirmed live.
Audited all70 legacy alignments against1,944 reference genes in the frozen
proteomes. Gene inventories/row lengths pass throughout;69 families preserve
ungapped residues exactly. RefOG023 matches only after removing23/11 stop
symbols from ENSP00000487059/ENSP00000486295. All34 reference stop symbols
occur in those two records;102 X residues also occur in the reference inputs.

Prepared a separate pinned MAFFT7.525 rebuild, leaving legacy outputs intact.
Explicitly remove stops, preserve X, force amino-acid mode, single thread per
family/eight concurrent families. Snapshot the entry script and libexec
companion files before/after execution; set the binary directory explicitly.
Validate normalized sequence preservation and calculate canonical-only mean
pairwise identity, with missing family values for any incomparable pair.
Per-family commands/status/failures are retained; no accuracy scoring occurs.
Fifteen focused tests passed. Full-suite and frozen submission details follow.

Full unit suite **1015 passed in25.48s** before freezing the alignment executor.

Alignment executor frozen**a753d0d** at
`benchmarks/work/publication_ob_reference_alignments_v1`; submitted**21309**
(eight CPUs32GiB/two hours). Output `benchmarks/results/ob_reference_alignments_v1`;
log `benchmarks/work/ob_reference_alignments_21309.log`. QfO21307 remains live.
Next terminal/native alignment admission and sequence-feature join, followed
by the frozen84-endpoint stratified outcome analysis. No family-error result
is claimed from this preparatory work; the full publication goal remains open.

### Native alignment admission and frozen family strata (2026-09-16)

Previous turn made progress: legacy alignment audit, tested normalized
rebuild and frozen21309 submission. Reread full objective;21307 remains live.
Alignment21309 completed0:0 in4:24; all70 families/1,944 distinct reference
genes succeeded. New assembler independently checked frozen executor/tool
inventories, scheduler completion, status/command/input/output provenance,
normalized sequence preservation and canonical-only identities. Verified the
entire251,378-row feature table against gene/species inventory, recomputing
all reference proteins' features from the frozen FASTAs.

Prepared outcome-independent family categories at
`benchmarks/results/ob_error_strata_v1/manifest.json`; preserved snapshots
`ob_reference_alignments_prepared_20260916.json` and
`ob_error_strata_prepared_20260916.json`. Size bins34/30/6; single/multi-copy
6/64; lower/higher/missing identity29/29/12; short-relative/not-short/missing
40/30/0; concentrated/not-concentrated/missing composition1/69/0.
Identity median among58 evaluable families:0.5783660968575601. Twelve missing
families have at least one pair without canonical overlap; do not assign zero
identity or silently omit them. The one-family composition bin receives no
bootstrap interval under the frozen minimum-five rule.

Illustrations were selected by SHA256 ranking of canonical RefOG filenames
(including .txt) within each occupied bin, with deduplication:005,014,021,024,
038,067. No method accuracy/error outcome entered selection. All70 families
remain scheduled for mechanistic tracing. Sixteen focused tests passed;
full unit suite **1031 passed in25.92s**. Next the84-endpoint stratified outcome
assembly, with fresh scoring/provenance crosschecks. Parameter/scaling/domain/
biological-application/QfO/release requirements remain open.

### Frozen 84-endpoint stratified scoring prepared (2026-09-16)

Previous turn made progress: full native alignment admission and frozen
family strata. Reread full objective;21307 is still live, no automatic restart.
Implemented a scoring assembler with fresh feature/alignment admission before
method outcomes, exact native-format parsers, input-universe checks, full
reference per-family sufficient-statistic reproduction, and official-score
crosschecks. Retains native prediction coverage without silently adding
unassigned genes. In retained OrthoFinder3.1.5 full output, Orthogroups.txt
is root-HOG-postprocessed (documented in the manuscript), not the raw MCL
checkpoint solely because of its filename.

Extended the shared paired bootstrap with an optional global multiplicity
count. Default results remain unchanged; the new panel retains84 endpoints
across all14 bins/two contrasts/three metrics. Empty bins have null estimates;
bins below five families have point estimates but no intervals. Weighted
counts and family wins/ties/losses are retained, with full-reference scoring
conventions unchanged within each bin. Sparse adjusted percentile tails
(about six draws each at20,000 replicates) are explicitly caveated.

Twenty-six focused bootstrap/assembly tests passed, covering global-adjustment
invariants, missing/small bins, malformed panels, native parsing, changed
inputs, mismatched sufficient statistics and official-score disagreement.
Full unit suite and frozen scoring submission details follow.

Full unit suite **1047 passed in26.90s**; focused assembler tests rerun after
adding job/scorer provenance and an explicit unknown-method guard.

Frozen scoring executora336255 at
`benchmarks/work/publication_ob_stratified_errors_v1`, job21310 completed0:0
in18s. All native predictions cover251,378 genes; fresh full-reference family
counts and official-score crosschecks exactly reproduce retained benchmarks.
All84 endpoint records regenerated from sufficient statistics; input/source
hashes rechecked. Snapshot `ob_stratified_error_results_20260916.json` SHA256
a8a072e9a3917fa7dabf25dba4affb8c3fd6b30c02286dda3efe882baceb3cac.
Eleven bins have66 interval endpoints, one bin is descriptive-only, two empty.
None of22 adjusted F1 intervals exclude zero. Two precision advantages/seven recall
deficits exclude zero across overlapping strata, not independent replication.
Both modes have lower recall in short-relative families; no fragment mechanism
is inferred. Full84-effect Markdown and manuscript/claim updates preserve
negative, neutral, missing and small-bin results.

QfO21307 completed0:0 in32:26. All three workers reproduce349,898 groups with
SHA2568c162782bba955078a8f7e4e6088df5b34c9f15b055803224e513fdea8b33cdd,
matching diagnostic21295 but not capture21305. All recorded modules,71 loaded
libraries, Python/packages/environment/host/platform/affinity agree; independent
file/hash/full-partition checks passed. Preserved `qfo_saved_graph_repeats_20260916.json`
SHA25633c3e27f3ec750335c6368de55b2240deb26f74b0864b186fab79761c13a1e75.
This is single-CPU configuration repeatability, not an identified root cause.
Next controlled one-versus32-CPU affinity repeats on the same graph/worker,
then process-context tests as needed; no default change or full-profile restart
based on selecting a preferred partition. Error-analysis figures/tracing and
the wider parameter/scaling/annotation/application/QfO/release tasks remain open.

### Controlled QfO CPU-affinity panel prepared (2026-09-16)

The preceding user-facing goal-prompt turn was no analysis progress. Reread the
full objective and revalidated current repository and scheduler state: no active
OrthoHMM jobs; unrelated job20915 remains untouched. Proceeding with the next
available reproducibility experiment, not shrinking the publication objective.

Added an optional four-arm affinity panel to the saved-graph observer while
preserving its default three-repeat interface and the earlier frozen executor.
Prespecified one/32/one/32 CPUs within one allocation, same graph/worker and
one-thread BLAS/OpenMP settings. Worker affinity is set before native imports;
inherited/requested/actual affinity and loaded software identities are checked.
Every partition is retained, with same-affinity and between-arm comparisons.
No accuracy scoring, preferred-partition selection, or scientific default change.
Protocol: `QFO_AFFINITY_DIAGNOSTIC_PROTOCOL_20260916.md`.

Nineteen focused tests passed, including real native-worker execution both with
and without explicit affinity, isolate preservation, parent affinity preservation,
allocation guards, and actual-affinity mismatch rejection. Full unit validation
and frozen submission details follow. Existing unrelated sample changes remain
untouched; their whitespace warnings are not part of this milestone.

Full unit suite **1055 passed in29.11s**. Scoped whitespace validation passed.

Committed/pushed executor **cc3bffa**, frozen at
`benchmarks/work/publication_qfo_saved_graph_affinity_v1`. Submitted job**21311**
from its frozen batch script,32CPU/64GiB/two-hour limit. Scheduler independently
confirmed RUNNING with zero restarts; no results yet. Output destination:
`benchmarks/results/qfo_saved_graph_affinity_v1`. Next verify all four native
worker records/partitions and interpret the controlled affinity contrast before
any full QfO replay. While it runs, error-analysis figures and mechanistic tracing
can proceed. All broader publication requirements remain active.

Push again reports21 dependency alerts (one critical, seven high, eleven moderate,
two low). These remain untriaged release work; frozen scientific environments
were not modified to address them during this diagnostic.

### Complete stratified error figure (2026-09-16)

Previous turn made progress: tested/pushed affinity executor and confirmed21311
live. Reread the objective and polled21311 again: still RUNNING, no restart.
While it runs, completed the OrthoBench error figure from the hash-pinned
a8a072e9 result snapshot. All14 strata/84 endpoints remain represented:
66 interval-bearing effects, six descriptive points,12 explicitly nonestimable
endpoints. F1/precision/recall share a horizontal scale; family counts, missing
rows, exploratory scope, overlapping strata and sparse adjusted tails remain
visible. No outcome-based filtering or biological-mechanism claim was added.

New `plot_ob_stratified_errors.py` validates feature-based membership, complete
methods/metrics, uncertainty specification, full-reference sufficient-statistic
point estimates, contrasts, interval nesting and small/empty-bin conventions.
Eleven focused tests passed including altered-input rejection and complete
rendered point/interval/missing-row inventory with text-bounds checks. Generated
PNG/PDF/SVG plus provenance manifest under
`figures_ob_stratified_errors_20260916`; PNG visually inspected with no clipped
labels or overlapping content. Integrated figure/caption into manuscript.
All70-family mechanistic tracing, independent domain/fragment/duplication
annotations, QfO strata and the wider publication requirements remain open.

Full unit suite **1066 passed in28.77s**. Independently rechecked result, plotter
and all three generated figure checksums against the manifest. Job21311 remains
RUNNING at5:36; first worker recorded actual affinity[8]. No terminal output or
affinity-effect conclusion yet. Next: complete mechanistic stage tracing and
admit the four-worker QfO result when terminal.

### All-family retained-stage tracing prepared (2026-09-16)

Previous turn made progress: complete84-endpoint figure and manuscript link
committed/pushed. Reread objective;21311 remains scheduler-confirmed RUNNING
at12:19, with no completed worker comparison yet. Preserved its live process.

Inspected retained checkpoints and actual merge-sidecar schema (a list of8440
events, not the older analysis script's expected wrapper object). Added a new
tracer for all70 families through six fixed checkpoints with full pair-level
directional cache evidence and group membership. Reconstructs candidate groups
from all logged merges, verifies complete input partitions and root-HOG candidate
boundaries, and crosschecks fresh full-reference scores at four frozen factorial
stages. Keeps raw co-membership descriptors distinct from official scoring and
records missing prefilter/edge/tree-level causal evidence explicitly.

Fifteen focused tests passed, covering native membership integrity, reference
versus unlabelled incident pairs, directional search values, nonmonotone grouping
transitions, invalid scores, iteration-start merge snapshots, and malformed or
nonreconstructing merge traces. Execution protocol extended without new endpoint
selection or inference changes. Full unit validation and frozen submission follow.

Full unit suite **1081 passed in29.42s**. Large pair-level output stays outside
Git; the committed result will retain its checksum and full70-family summaries.

Frozen3cb213c executor at`publication_ob_family_trace_v1`, job21312,
FAILED1:0 in1s before output creation: an incorrect new disjoint-reference
assumption. Independent input audit confirms1,945 memberships/1,944 unique
genes: FBpp0309618 belongs toRefOG021 andRefOG068. Existing scorers preserve
these references; no previous score change is indicated. Removed the tracer's
unsupported disjointness restriction and added an explicit overlap inventory
and regression test. Both assignments and all70 families remain intact.
Preserved failed job/log/frozen source; corrected execution uses separatev2
paths and is not an automatic retry. Full validation and job ID follow.

Corrected full unit suite **1082 passed in28.09s**, scoped whitespace clean.
QfO21311 first one-CPU arm completed with the prior single-CPU partition SHA256
8c162782bba955078a8f7e4e6088df5b34c9f15b055803224e513fdea8b33cdd;
remaining arms still required before interpreting the affinity contrast.

Corrected executorcdea3c3 committed/pushed and frozen at
`benchmarks/work/publication_ob_family_trace_v2`; job**21313** submitted with
oneCPU/32GiB/one-hour limit and no automatic requeue. Scheduler confirms
RUNNING at11s, empty stderr so far. Output target`benchmarks/results/ob_family_trace_v2`.
Next admit terminal results against all source/pair-table hashes, independently
verify family transitions/merge reconstruction, and inspect the frozen illustrative
families. Both21313 and21311 remain live; no publication-completion claim.

### Native pair-trace admission and branch correction (2026-09-16)

Previous turn made progress: corrected overlap handling and frozen21313
submission. Reread objective;21313 is now COMPLETED0:0 in41s,21311 still live.
All70 families,40,733 pair rows,1,945 memberships/1,944 distinct genes retained.
All8,440 candidate merges reconstruct the candidate partition;168 logged events
touch reference genes. Pair table4,390,860bytes remains outside Git, SHA256
aabf71d4f79aff06b18890b49996895cd729bd6ac1f6fb0b7f14e79ab5295e4b.

Frozen replay source review found a reporting error in the new tracer: checkpoints
were treated as a linear chain although profile expansion starts from unrefined
multipass clusters. Preserved original82883d3d result; corrected source-defined
branches and separately labeled refined-endpoint comparison. No pair membership,
inference output, score or method default changed. Added a native-file admission
that verifies exhaustive pair coverage, exact booleans/species/membership,
fresh full-reference scores and all group summaries/source hashes. Tests reject
missing/duplicate/reversed/malformed pairs and wrong native memberships.

Admission output`ob_family_trace_admitted_v2` adds helper provenance; first
admission output is preserved. Committed snapshot
`ob_family_trace_verified_20260916.json` SHA256
bda00fb593b357bc8f07e43544feae598150ae891ed82c988fb339a92349de0c.
All-family table, full six frozen illustrations, interpretation and manuscript
text retain adverse and neutral cases. Profile-on versus off refined endpoints
gain836/lose160 raw within-family pairs; candidates gain4,272; final rootHOGs
lose1,575. These are descriptive co-membership counts, not official recall or
causal attribution. Rejected-hit/added-edge/gene-tree/constraint mechanisms,
independent annotations and biological application remain open.

Focused26 tests passed; full unit suite **1092 passed in29.51s**.21311 remains
RUNNING at21:36, first one-CPU arm only completed; no affinity conclusion yet.

### Recorded root-lineage versus constraint decomposition (2026-09-16)

Previous turn made progress: admitted complete native pair trace and corrected
branch semantics. Reread objective;21311 still running. Source inspection
established the node table retains pre-constraint reconciliation calls, while
rootHOGs are post-constraint. Frozen phylogeny.py/pipeline.py are byte-identical
to reviewed main sources. Added a reconstruction helper plus full native-file
assembler, supporting only the frozen species_overlap/positive_paralogy rules.

All250 reference-incident candidate families reconstructed exactly:114 have
node tables,136 valid bypass cases. All193 internal constraints included,
121supported/72detached. This scope differs from168 directly reference-incident
logged merges and does not claim the genome-wide8,440-constraint audit.
Node topology, child membership/species overlap, pair-event rules, selected
tree/checkpoint hashes, final candidate boundaries and all final group memberships
checked. Complete pair decomposition matches the prior40,733-row trace.

Native pair dispositions:16,264 already in different candidates;134 root-lineage
splits;1,441 subsequent constraint splits;22,894 retained. Four reference families
have root-lineage losses and six constraint losses, withRefOG011 shared. All61
other families have no post-candidate within-reference loss. Frozen illustrations
014and021 lose95/806 pairs at the constraint step; four other illustrations are
neutral at these steps. No biological truth or default change follows; the
unconstrained control's precision cost remains relevant.

Preserved exploratory extracts v1/v2; their only candidate-report difference
was unordered final-group serialization. Canonicalized group order forv3;
all candidate partitions and reference classifications agree after normalization.
Final564 provenance records independently rechecked. Snapshot
`ob_reconciliation_trace_20260916.json` SHA256
4d5822db16b7c326ce2ac923ecd98395ea4bf9d0ef77cfa0a3cae4a9d1d0ebb9.
Result interpretation, all affected families, six frozen illustrations, protocol
and manuscript updated. Fifteen focused tests include native-reconciler parity;
full unit suite **1107 passed in29.41s**.

QfO21311 remains RUNNING at34:10. First one-CPU and32-CPU arms both yielded
8c162782bba955078a8f7e4e6088df5b34c9f15b055803224e513fdea8b33cdd,
with recorded software identity excluding affinity equal. Two repeats remain;
do not infer affinity causation or general determinism from this partial panel.
Tree-history validation, edge/rejection tracing, independent annotations,
QfO/control/parameter/scaling/application/release requirements remain active.

### Prespecified parameter neighborhood prepared (2026-09-16)

Previous turn made progress: recorded reconciliation/constraint decomposition
tested and pushed. Reread objective and moved to the outstanding limited
parameter-robustness requirement. Fixed six one-at-a-time20% variants:
CPM0.08/0.12 around0.1; candidate min_norm0.024/0.036 around0.03;
candidate min_margin1.2/1.8 around1.5. Unchanged control required. This is an
exposed-data sensitivity panel, not adaptive optimization or new-default selection.
Protocol`PARAMETER_NEIGHBORHOOD_PROTOCOL_20260916.md` specifies inferred-tree
full-pipeline endpoints and18 planned OrthoBench contrasts/metrics, no denominator
reduction after failures. QfO extension requires a reproducible baseline first.

Implemented candidate-only preparation for the control and four threshold arms,
with fixed byte-admitted HMM/profile seed groups and cached hits. Scoped
analysis-only call overrides leave frozen production code unchanged, verify
exactly one engine call, restore the engine on failure, and retain nominal
wrapper reports separately from effective parameters. Unchanged control must
match both full candidate and merge-trace bytes before variants proceed.
No reconciliation or accuracy is evaluated by this preparer; CPM arms remain
separate required work. Ten focused override/restoration/guard tests passed;
full unit validation and submission follow.

QfO21311 still RUNNING at41:50. Third arm(one_cpu_1) produced partition
def06c421941743617dd3540a14f8d9dc4f19eeedd34c1eb0c86c92160500488,
unlike the first one-CPU arm8c162782. Recorded software identity and affinity
match between these one-CPU workers. Thus the interim panel no longer supports
affinity alone as an explanation or general single-CPU repeatability. Final
arm/full admission pending; no selection, restart or accuracy scoring.

Full unit suite **1117 passed in27.74s**. Scientific defaults and held-out
evaluation remain unchanged.

Frozen candidate preparer **a09b740** committed/pushed at
`benchmarks/work/publication_ob_candidate_neighborhood_v1`, job**21314**
submitted oneCPU/32GiB/one-hour limit/no-requeue. Scheduler confirmed RUNNING
at14s; output`benchmarks/results/ob_candidate_neighborhood_v1`. No accuracy
evaluation or completed downstream parameter panel claimed.

QfO21311 subsequently reached COMPLETED0:0 in43:09. Final manifest status
affinity_panel_complete: one_cpu_0/all_cpus_0/all_cpus_1 all8c162782;one_cpu_1
def06c42. Recorded software identity excluding affinity matches throughout;
both single-CPU workers also match affinity. Complete independent native-file/
partition admission is next, before committing a result snapshot or promoting
any reproducibility claim. This bounded panel does not isolate a root cause;
same-affinity disagreement remains unexplained and preserved without selection.

Before closeout21314 completed0:0 in1:12; manifest reports all five candidate
arms prepared unscored. Unchanged control passed both byte-equivalence gates.
Reported candidate counts:control54,445;norm_low54,370;norm_high54,540;
margin_low52,912;margin_high55,434. Independent admission and full-pipeline
reconciliation remain next; these counts are not accuracy results.

### Independent QfO affinity-panel admission (2026-09-16)

Previous turn made progress: frozen parameter panel/preparation executed and
both jobs reached terminal success. Reread objective and independently admitted
completed21311. New validator checks the exact four-arm plan, native snapshots,
execution records, same graph, actual/inherited/requested affinity, module/library
identities and fixed worker parameters. Recomputed all13 stored full partition
comparisons and checked248 unique provenance files. No failed check or partial
panel was treated as complete.

Verified same-CPU difference:349,898 versus349,950 groups;4,427/4,479 unmatched
groups, all976,504 genes present. Both32-CPU repeats are byte-identical349,898.
Same-affinity disagreement is now directly verified, not merely a summary-file
observation. Preserved snapshot`qfo_affinity_verified_20260916.json` SHA256
a927b00477f12e0a4c6548271942bc2f78d047bf0f7d0304d8aa0abd11954ec0.
Diagnostic and manuscript corrected to avoid a general single-CPU repeatability
claim. CPU availability alone is insufficient; effective native graph/state
still needs investigation. Static installed-binding inspection confirms the
Python wrapper calls optimizer.set_rng_seed with the supplied seed, but does
not identify a native cause.

Ten focused admission tests passed; full unit suite **1127 passed in26.82s**.
Next record actual igraph edges/weights and effective optimizer arguments in
bounded fresh workers; preserve every partition, with no accuracy selection.
Candidate-neighborhood21314 independent admission/downstream inferred-tree
execution, CPM variants and broader publication requirements remain active.

### Native optimizer-boundary probe preparation (2026-09-16)

Previous response supplied a goal prompt but did not advance experimental state;
this continuation reread the objective and resumed the available diagnostic work.
Added optional observation of the actual igraph object passed to find_partition:
ordered canonical undirected endpoints, ordered float64 weights, vertex/edge
counts, partition class and bound effective arguments including library defaults.
Chunked hashing preserves edge order, parallel edges, loops and isolates. Compare
against the saved arrays before optimization and recheck graph contents afterward.
The original call receives unchanged arguments. No scientific core edits.

Prespecified three fresh one-CPU workers using the same saved QfO graph, CPM0.1,
seed4 and one-thread native settings; all partitions retained, no accuracy scoring.
Instrumentation can change native memory/runtime context, so this is a boundary
diagnostic, not proof about uninstrumented execution or controlled timing.
Batch qfo_native_boundary_batch_20260916.sh targets a frozen executor worktree.
Full unit suite1133 passed25.67s before adding two additional real-worker probe
combinations; updated worker suite21 passed3.10s. Six boundary tests cover native
defaults, graph equality, ordering/weights, failures and restoration. Actual QfO
execution and independent result admission remain pending at this milestone.

Executor b277e5b committed/pushed and frozen at
benchmarks/work/publication_qfo_native_boundary_v1. Submitted no-requeue Slurm
job21315 (one CPU,64GiB,two-hour limit); scheduler confirmed RUNNING. Output is
benchmarks/results/qfo_native_boundary_v1. Next inspect terminal accounting,
independently compare native boundary records and all partitions, then use the
evidence to narrow graph-conversion versus optimizer-state explanations.

### Candidate neighborhood independently admitted (2026-09-16)

Previous turn made progress: native-boundary probe implemented/tested/pushed and
job21315 launched; scheduler reconfirmed RUNNING during this continuation.
Independently admitted completed21314 while that diagnostic runs. Exact five-arm
plan, nominal/applied parameter separation, one engine call, all70 provenance
records, full251,378-gene coverage, control byte equivalence and complete merge
reconstruction passed. Candidate counts control/norm_low/norm_high/margin_low/
margin_high:54,445/54,370/54,540/52,912/55,434; corresponding reconstructed merges
8,440/8,515/8,345/9,973/7,451 from62,885 fixed HMM seed groups.

Snapshot ob_candidate_neighborhood_verified_20260916.json SHA256
38a4cc8c16e8e45800c88af618f9b5a9b3d77ab6b3d0e5b79fa555b88e849514.
Admission verifies retained memberships and provenance, not independent search
or threshold-decision recomputation. No reference scores evaluated. Twelve new
tests; full unit suite1,147 passed27.64s. Four downstream inferred-tree variant
runs and two CPM variants remain pending, followed by prespecified18-endpoint
scoring. QfO21315 remains running; no native-boundary conclusion yet.

### Candidate neighborhood inferred-tree executor (2026-09-16)

Previous turn made progress by independently admitting the five candidate arms.
This continuation reread the objective and confirmed21315 still RUNNING. Added
four-arm downstream executor using admitted snapshot38a4cc8c, frozen factorial
scientific launcher and environment. Changes only candidate/constraint paths and
destinations, plus baseline checkpoint-source reuse. Keeps species-tree-mode infer
and all other frozen settings. Baseline native artifacts checked before/after;
all preparation provenance and environment rechecked. No reference scoring.

Native checkpoint code checks exact gene membership, sequence/config input hashes
and raw-tree checksums; species-tree cache checks selected marker inputs/config
and tree checksum. Therefore a genuinely unchanged inference input may reuse its
tree, but no baseline tree is supplied to a changed marker input. Reused source
artifacts remain protected by native atomic replacement and postflight identity
checks. Failures preserved without retries; successful execution still requires
independent native admission before scoring.

Batch defines four tasks,32CPUs/64GiB/four-hour limit each,max two concurrently.
Shared-machine incremental times are not controlled end-to-end scaling evidence.
Eleven focused executor tests passed, bash syntax passed, full unit suite1,158
passed29.27s. Freeze/push executor before submitting; CPM variants and complete
six-variant18-endpoint analysis remain pending.

Executor e8e86d8 pushed and frozen at publication_ob_candidate_phylogeny_v1.
Submitted no-requeue array21316: tasks0(norm_low)/1(norm_high) confirmed RUNNING;
tasks2(margin_low)/3(margin_high) PENDING for JobArrayTaskLimit. Output root
benchmarks/results/ob_candidate_neighborhood_phylogeny_v1. QfO21315 remains
RUNNING. Next validate terminal native outputs/coverage/checkpoint semantics,
prepare both CPM variants, and complete prespecified scoring without selection.

### CPM neighborhood replay executor (2026-09-16)

Previous turn made progress: frozen candidate-phylogeny array21316 submitted.
Reread objective;21316_0/1 and QfO21315 confirmed RUNNING,21316_2/3 queued.
Prepared control0.1 followed by prespecified CPM0.08/0.12, sequentially in one
32CPU/64GiB/four-hour job. Uses the exact hash-verified successful OrthoBench
replay command/source/runtime from21138, changing only resolution/destinations.
All four unchanged-control stage partitions must match baseline bytes before
either variant starts. RBNH grouping, singleton assignment, profile construction,
profile search/reclustering and refinement rebuilt from fixed cached search hits.
No reused baseline profile-expanded groups and no reference scores.

Output captures stage partitions, profile metrics, commands, input/core/runtime
provenance and per-arm time logs. Candidate expansion and inferred phylogeny for
these two CPM variants remain subsequent requirements; do not call replay alone
the complete robustness experiment. Shared-node incremental timing excludes the
original all-to-all search. Thirteen focused command/gate tests passed, batch
syntax passed, full unit suite1,171 passed30.31s. Freeze/push before submission.

QfO21315 first instrumented optimizer call has returned with recorded976,504
vertices/24,148,515 edges and expectedCPM0.1/seed4/default2 iterations. This is a
partial live observation only; full repeat comparison and independent admission
still pending. No determinism or root-cause claim from one worker.

Executor1bae2d2 pushed/frozen at publication_ob_cpm_neighborhood_v1; no-requeue
job21319 submitted and confirmed RUNNING. Output benchmarks/results/
ob_cpm_neighborhood_v1. Existing21316_0/1 and21315 remain RUNNING;21316_2/3
remain queued under the array concurrency limit. No job restarted or superseded.

### QfO pre-optimizer endpoint mismatch (2026-09-16)

Previous turn made progress: CPM executor/tests/push/submission21319. Objective
reread. Scheduler now confirms21315 FAILED1:0 in11:49. First worker completed;
second failed the explicit native graph equivalence gate before optimization.
Same edge/vertex counts and ordered weights, different ordered endpoint hash:
native cb8777a1... versus saved dd9c0c07... . Third worker never launched. Full
failure and partial results preserved, not admitted as a successful repeat panel.
Fresh saved-array fingerprint and four input file hashes verified; recorded
software/environment/affinity identical between workers. See updated drift
diagnosis and two committed partial-evidence snapshots. This narrows this run's
discrepancy to before optimization but does not prove which conversion step or
whether instrumentation caused it, nor explain all historical variation.

Enhanced optional probe records differing-edge counts/examples and compares the
actual frozen caller's graph_edges array against both native and saved endpoints
when a mismatch occurs. Observer only; no scientific core changes. Separate v2
batch/output to investigate mechanism without replacing failed21315. Two added
tests cover constructor/native separation and caller-array capture;29 focused
probe/worker tests passed. Full unit suite1,173 passed29.84s. Candidate array21316
and CPM21319 remain running; no scores or default changes.

Enhanced observer executor984355c pushed/frozen at publication_qfo_native_boundary_v2.
Submitted no-requeue job21321, separate output qfo_native_boundary_v2; scheduler
confirmed RUNNING. Original21315 remains terminal failed with all artifacts.

### Normalized-support phylogeny variants admitted (2026-09-16)

Previous turn made progress by preserving native mismatch evidence and launching
the detailed QfO probe21321. Objective reread; both normalized-support tasks of
21316 completed0:0 (low7:35/high6:05), margin variants still running. Added
independent native validator with exact array/raw-job identity, frozen command
and executor, preparation/environment hashes, postflight, native metadata,
constraint accounting, input/outputs/tree provenance, species coverage and
complete nonoverlapping root-HOG coverage checks. Execution artifact hashes
recomputed:63,779low/64,627high. Fresh baseline audit confirms scientific source
artifacts unchanged. Validator adapts only candidate/constraint records for the
shared native checker, preserving all other baseline requirements.

Initial validation attempt rejected whole baseline admission equality because
the verifier path differs between frozen executor and main checkout. Corrected
to compare native manifests, metrics, species tree, partition, membership and
status instead of observer location; regression test added. No native outputs
changed or inference rerun. This validation attempt produced no admitted file.

Both admitted variants retain251,378 genes, with zero cross-candidate merges:
norm_low54,370candidates/59,729rootHOGs;8,515constraints=5,903supported+2,612detached.
norm_high54,540candidates/59,822rootHOGs;8,345constraints=5,782supported+2,563detached.
These are integrity/count results, NOTaccuracy or robustness conclusions.
Snapshots ob_candidate_norm_low_native_verified_20260916.json SHA256
2f9c02375fcf7a40c3d11237193dc329545442e199cbb2140df714a4e5bcc6f9;
ob_candidate_norm_high_native_verified_20260916.json SHA256
3d81a3bf0c337f8703104abe2260651be29b24af54b6795cb57db3f48fdb3a51.

Fourteen focused validator tests passed; full unit suite1,187 passed34.43s.
Manuscript limitation updated with partial QfO pre-optimizer discrepancy, without
claiming a proven cause. Margin-variant admission, CPM replay/candidate/phylogeny,
six-variant scoring and QfO21321 diagnostic conclusions remain pending.

### Fixed parameter-panel statistics and control replay (2026-09-16)

Previous turn made progress by admitting both normalized-support phylogeny
variants. Objective reread, remaining jobs confirmed RUNNING. Added the fixed
statistical component for the six-variant parameter panel: exact complete arm
accounting, successful baseline required,20,000 paired RefOG draws/seed20260918,
18 planned endpoints even with terminal failures. Missing/extra/duplicate arms,
undocumented failures and mismatched reference families rejected. Failed variants
receive explicit missing intervals, never zero accuracy. If all variants fail,
retain baseline statistics and failure reasons without fabricated bootstrap CIs.
Native admission and official-score verification remain caller requirements;
this helper does not replace those gates or constitute the final assembler.
Nine focused tests; full unit suite1,196 passed33.64s. No panel scores calculated.

CPM21319 unchanged control completed its four stages and passed byte/partition
equivalence (51,181multipass;63,245refined;50,894profiles;62,885profile-refined).
All four retained stage hashes independently rechecked after gate observation.
CPM0.08 is now running;CPM0.12 remains planned after it. This control success is
not completion of the whole CPM panel or evidence of universal repeatability.
Margin-variant21316_2/3 and detailed QfO21321 remain running. Next finish/admit
these outputs and CPM downstream candidate/phylogeny before six-variant scoring.

### Sequence-search figure and high-margin admission (2026-09-16)

Previous turn made progress with fixed statistics implementation and verified
CPM control. Objective reread and remaining jobs checked. Generated missing
sequence-search control figure from hash-pinned b3740104... admitted results.
Three-arm score table plus all six paired F1/precision/recall effects, same x-axis,
nominal and multiplicity-adjusted intervals. Clearly labels disabled profile/
candidate/phylogeny stages, unmatched sensitivity/calibration and no established
HMM F1 advantage. Plot validation recomputes scores from sufficient statistics,
checks coverage, contrast arithmetic, interval nesting and exact endpoint plan.
PNG visually inspected: nonblank, legible, no clipping or overlapping labels;
PNG/PDF/SVG plus provenance manifest retained and manuscript linked.
Ten focused tests passed; full unit suite1,206 passed29.20s.

Completed21316_3 (margin_high)0:0 in11:37 independently admitted via the existing
native validator;75,035 execution artifacts verified. All251,378genes retained,
55,434candidates/60,156rootHOGs, zero cross-source merges. Constraints7,451=
5,361supported+2,090detached. Unscored snapshot
ob_candidate_margin_high_native_verified_20260916.json SHA256
07a0035414254da6ccbb5fbc8928fd22720caedf48e77188ba615d220cc1e4ec.
Only margin_low remains running in21316. CPM21319 and QfO21321 still running;
first v2 QfO worker matches8c162782... partition, not a completed repeat panel.
No parameter accuracy results or new defaults inferred from these observations.

### All candidate-threshold native variants admitted (2026-09-16)

Previous turn made progress with sequence-control figure and high-margin native
admission. Objective reread. Last candidate task21316_2 completed0:0 in19:33;
fresh native validation passed75,065 execution artifacts, all251,378genes,
52,912candidates/59,329rootHOGs and zero cross-candidate merges. Constraints9,973=
6,527supported+3,446detached. Snapshot
ob_candidate_margin_low_native_verified_20260916.json SHA256
94cfce4a7b101a561eb023c5c0cec51f76add0be7f841fad0d282081aaa762e2.
All four threshold variants now independently admitted, still unscored.

Added complete-panel CPM replay admission code for21319: require terminal success,
exact three-arm plan/executor/commands/input records, frozen source/runtime,
complete stage partitions, native profile-build/iteration accounting and fresh
four-stage control equivalence. Positive profile-build evidence is mandatory for
this panel; it is not independent validation of every profile score. Eleven
focused tests passed; full unit suite1,217 passed28.83s. CPM0.08 metadata/profile
accounting checked on real completed output; full panel admission deliberately
awaits still-running CPM0.12. Candidate expansion and inferred phylogeny for both
CPM variants remain pending. Detailed QfO21321 still running, with no new completed
panel or root-cause conclusion. No default promotion or accuracy selection.

### CPM replay admitted and candidate preparation ready (2026-09-16)

Previous turn made progress by admitting all candidate-threshold phylogeny arms
and implementing CPM admission. Objective reread;21319 COMPLETED0:0 in24:45.
Independent full-panel admission passed77 provenance records, all12 stage
partitions, control equivalence, exact native command/parameters/runtime and HMM
profile-build accounting. Refined seed groups control62,885/CPM0.08 62,895/
CPM0.12 62,733. Snapshot ob_cpm_replay_verified_20260916.json SHA256
ceb8f9fe7c12b35317cde99c3d027d4d1080d398403afc04e57fb029ce1b460a.

Prepared frozen-core candidate executor using each arm's own admitted HMM-refined
seed with unchanged satellite_v2 parameters; same cached hits, complete coverage,
merge reconstruction and control candidate/trace byte-equivalence gates. Six
focused tests, bash syntax, full unit suite1,223 passed27.64s. No reference scores.
Next freeze/push/submit candidate preparation, independently admit it, then run
the two CPM inferred-tree variants before complete six-variant scoring.

During this continuation21321 terminatedFAILED1:0 in22:26. First two workers
completed; third stopped before optimization on a native endpoint mismatch.
Detailed observer reports constructor int32 C-contiguous array identical to saved
endpoints, but six native edges differ at indices23,493,880..23,493,885. This
narrows observed divergence past the Python constructor-input array to native
graph construction/storage/access, not yet a proven library defect or explanation
of every historical partition. Preserved qfo_native_constructor_mismatch_20260916.json
SHA25650bd555569203a477567f8bbca1e29718b82c42746c62bace953f9e39764c5c6;
full failed run remains benchmarks/results/qfo_native_boundary_v2. No retries or
default fixes yet. Independent witness checks and constructor-only experiments next.

CPM candidate executor070051d pushed/frozen at publication_ob_cpm_candidates_v1;
no-requeue job21322 submitted and confirmed RUNNING (one CPU,32GiB,one hour).
Output benchmarks/results/ob_cpm_candidates_v1; no phylogeny or scoring yet.

### Construction-only QfO diagnostic prepared (2026-09-16)

Previous turn made progress with CPM replay admission, candidate executor and
submission21322. Objective reread. Installed igraph Python wrapper inspected:
NumPy arrays are converted through numpy_to_contiguous_memoryview, which calls
numpy.require at native igraph integer width. This identifies a concrete input
conversion to test, not a proven defect. Added six-worker alternating diagnostic:
three original int32 inputs and three explicit int64 copies, identical frozen
graph handling, no optimizer invocation or partition/scoring. Records native
fingerprints, saved/original/converted array differences and bounded mismatch
witnesses via tuple/source/target/get_eid. Explicit copy changes allocation as
well as dtype; a clean panel cannot establish universal correctness.

First fresh-worker tests exposed a missed hook (function-local igraph import);
corrected to instrument Graph.__init__ while preserving the Graph class. Both
fresh-worker modes now prove observation occurs and no partition is generated.
Full unit suite1,226 passed31.12s; focused three tests passed again after adding
converted-array verification. Batch syntax passed. Scientific core unchanged.
Freeze/push before submitting bounded one-CPU64GiB/one-hour diagnostic.

CPM candidate21322 completed0:0 in1:23. Preparation reports control54,445groups/
8,440merges;CPM0.08 55,349/7,546;CPM0.12 53,548/9,185, each with full coverage
and merge reconstruction. These are unscored executor observations, not yet
independently admitted. Next admit them and run two inferred-tree variants.

Construction executor f2827a6 pushed/frozen at publication_qfo_construction_v1;
job21323 submitted no-requeue and confirmed RUNNING. Separate output
benchmarks/results/qfo_construction_v1. Previous failed boundary runs preserved.

### CPM candidates independently admitted; phylogeny ready (2026-09-16)

Previous turn made progress by launching bounded construction-only diagnostic.
Objective reread,21323 confirmed RUNNING. Independently admitted completed21322:
hash-pinned preparation/replay, own CPM HMM seeds, unchanged satellite parameters,
complete input coverage, logged merge reconstruction, count consistency and
control candidate/trace byte equivalence all pass. Candidate groups control54,445,
CPM0.08 55,349 andCPM0.12 53,548; all251,378genes. Snapshot
ob_cpm_candidates_verified_20260916.json SHA256
5acd1c56fe72e6267913a1170efca7f5f189960785422b7194f6a5fcc545a5fd.

Existing candidate phylogeny executor extended with explicit --cpm panel selection,
separately pinned admission and output directory. Commands retain inferred species
trees, native exact-input raw-tree reuse and frozen reconciliation/constraints.
No baseline threshold-panel behavior changed. Two32CPU/64GiB/four-hour tasks,
max two concurrent. Full unit suite1,236 passed30.04s;21 focused admission/runner
tests passed; batch syntax passed. No accuracy scores or default promotion.

Next freeze/push/submit both CPM phylogeny arms, then independent native admission
and complete six-variant scoring. QfO construction diagnostic first int32 worker
matched saved arrays; remaining workers running, no dtype conclusion yet.

Executor7e8c3e1 pushed/frozen at publication_ob_cpm_phylogeny_v1. Submitted
no-requeue array21324; tasks0(CPM0.08)/1(CPM0.12) both confirmed RUNNING.
Output benchmarks/results/ob_cpm_phylogeny_v1. QfO21323 remains RUNNING.

### CPM native admission prepared (2026-09-16)

Objective reread. Previous prompt-only turn did not advance authoritative project
state; resumed implementation after confirming21323 and both21324 tasks RUNNING.
Extended native candidate validator with explicit --cpm selection: pinned CPM
admission, executor7e8c3e1, array21324 and separate output directory. Threshold
panel retains its original executor/job defaults. Adapter retains each CPM arm's
own seed provenance rather than the baseline seed. Tests reject cross-panel
scheduler and execution identities;41 focused tests and full1,248 unit tests
passed (36.84s). Scoped diff check passed; unrelated sample-output whitespace
reported by repository-wide diff check left unchanged. Scientific core unchanged.

QfO construction21323 remains running. Four completed worker observations in its
progress report include an explicit_int64 worker with six native endpoint
differences at indices23493880..23493885. Both original int32 and converted int64
arrays match saved endpoints; native tuple/source/target agree on mismatches and
get_eid returns -1 for expected pairs. Explicit conversion alone therefore does
not eliminate the observed failure. This is partial diagnostic evidence, not an
admitted completed panel or proof of a particular library defect. No optimizer
or accuracy scoring invoked by this probe. Preserve all workers and await terminal
report before full admission. Next: admit CPM phylogeny after terminal success,
then assemble all six prespecified parameter contrasts with fixed multiplicity.

### Six-variant scoring assembler prepared (2026-09-16)

Previous turn made progress via tested/pushed CPM native admission support.
Objective reread;21323 and both21324 tasks remain scheduler-confirmed RUNNING.
Added assemble_ob_parameter_neighborhood.py: fresh baseline and all six native
admissions must succeed before reference access; pending/invalid runs abort the
complete-panel assembler. Exact arm identities, complete FASTA coverage, native
root-HOG provenance, official OrthoBench score agreement and post-score input
integrity are required. Existing statistics helper retains20,000 paired RefOG
bootstrap replicates/seed20260918/fixed18 endpoints; no best-arm selection.
JSON includes native admission, conversions, scoring helper hashes, coverage and
incremental phylogeny resource measurements. Reports explicitly exclude upstream
costs and do not interpret cached shared-node timings as end-to-end comparisons.

20 focused assembler/statistics tests passed; full1,259 unit tests passed38.53s.
After adding helper-hash tracking, all11 assembler tests passed again. No actual
parameter scores calculated while CPM native jobs remain pending. QfO progress
now contains five construction-only observations: four matching, one explicit
int64 worker with six native mismatches despite intact source/converted arrays.
Final worker/report still pending; no optimizer or accuracy evaluation in probe.
Next: admit terminal CPM runs and execute complete-panel scorer; independently
inspect final construction probe and continue diagnosing native endpoint failure.

### Construction evidence admitted; CPM low admitted (2026-09-16)

Previous turn progressed via the tested/pushed six-variant scoring assembler.
Objective reread. Construction-only21323 completed0:0 in18:44. All six workers
preserved: five match; explicit_int64 repeat1 has six native endpoint mismatches
despite intact original and converted arrays. New independent admission checks
terminal state, frozen executor/runtime, exact worker inventory/software/inputs,
257 file records, native-file/report agreement and complete reconstruction of
each ordered native endpoint hash from saved arrays plus bounded witnesses.
Snapshot qfo_construction_verified_20260916.json SHA256
499f90e06b6da311356c80f435149ba0cf132f75b13e2eacf720e02ea8f2382e.
Documented evidence and limitations in QFO_CONSTRUCTION_DIAGNOSTIC_20260916.md.
Explicit int64 conversion alone is not a sufficient fix; cause remains unresolved.
No optimizer execution, partition selection or accuracy evaluation in this probe.

Initial admission implementation had a syntax typo caught at import, corrected
before any output was written.14 focused tests then passed, including native
hash reconstruction and corrupt evidence rejection; full1,273 unit tests passed
29.44s. Native scientific code and frozen runtime unchanged.

CPM low21324_0 completed0:0 in16:42 and independently passed native admission:
55,349 candidate families,60,635 root HOGs,251,378 genes fully preserved,
zero cross-source merges,1,944 split families;7,546 constraints with5,178
supported and2,368 detached. Snapshot ob_cpm_low_native_verified_20260916.json
SHA2568f786c286894fc6dab757427048d3e1182b425a9dfaf1090dc71271bc63f071e.
CPM high21324_1 completed0:0 in17:44; native validation now running. No accuracy
scores yet. Next complete its admission and execute all six parameter contrasts.

CPM high admission now passed:53,548 candidate families,58,815 root HOGs,
251,378 genes preserved,zero cross-source merges,2,042 split families;
9,185 constraints with6,525 supported and2,660 detached. Snapshot
ob_cpm_high_native_verified_20260916.json SHA256
f6749ff370ebfb6f83318ec376f96db839958e0cd00da053e8d66b3449d21a16.
Both CPM and all four threshold variants are now natively admitted, still
unscored. Complete-panel accuracy/uncertainty assembly is the next action.

### Complete OrthoBench parameter panel scored (2026-09-16)

Previous turn made progress by admitting both CPM outputs and construction
evidence. Objective reread. Executed complete-panel assembler from main87d3b51;
all seven native outputs freshly admitted before references, all251,378 input
genes retained, and all scores independently crosschecked with official scorer.
No pending/failed variants.20,000 paired RefOG draws,seed20260918,fixed18
planned F1/precision/recall endpoints. Report SHA256
20adcd019644aa59e97a03cc951d66d900c9f3d705f79b67b798e52ae06d509b,
snapshot ob_parameter_neighborhood_results_20260916.json; generated full table
OB_PARAMETER_NEIGHBORHOOD_RESULTS_20260916.md.

F1 control74.106074;CPM low74.586147/high71.463468;norm low74.008570/high
74.110826;margin low74.973144/high73.274559. All six adjusted F1 intervals
include zero. Only CPM-high recall excludes zero after adjustment:delta-2.010088,
CI[-6.291062,-0.147598]. No new default or superiority claim. Norm changes
affect1/3families respectively, not evidence of general invariance. Family
wins/ties/losses and incremental cached shared-node resources retained in JSON.

Added parameter figure generator that recomputes all paired statistics from
sufficient counts before rendering. Eight focused tests passed; full1,281 unit
tests passed32.32s. Generated PNG/PDF/SVG plus manifest; visually inspected PNG:
all seven table rows,18 effects and caveats legible, no overlap/clipping. Common
axis preserves comparability; tiny norm effects remain numerically documented.
Manuscript updated with parameter results and independently admitted QfO
construction findings. QfO reproducibility, broader independent/biological
validation, matched scaling and portable release/archive work remain unfinished.

Staged diff check reports Matplotlib-generated SVG path whitespace only; keep
generated bytes unchanged to preserve the figure manifest. Source/document/JSON
checks excluding that generated SVG pass.

### Direct graph-stage diagnostic ready (2026-09-16)

Previous turn progressed by completing OrthoBench parameter scoring and figures.
Objective reread. Next unresolved primary-benchmark prerequisite is native QfO
graph reproducibility. Added probe_qfo_direct_graph.py: six fresh one-CPU workers
alternate minimal imports and frozen-worker imports, three repeats each, using
the same admitted saved RBNH arrays and original int32 constructor path. Direct
construction does not execute the scientific worker or Leiden optimizer. Inspect
all native endpoints immediately after construction and after weight assignment;
compare original constructor inputs to saved endpoints at each stage, preserve
bounded witnesses and full weighted graph fingerprint. Record imports, libraries,
source/input hashes, environment, affinity and partial failures. No accuracy
scores, clustering partitions, retries or new defaults.

This intentionally omits frozen worker bookkeeping; instrumenting the pre-weight
stage changes allocation and timing. Import-mode differences are not equivalent
to a matched historical replay, and matching observations cannot prove general
correctness. The aim is to localize an observed failure before optimization, not
choose a favorable partition. Previous construction panel and all failures stay
unchanged. New input gate pins its admitted snapshot499f90e0... and rechecks all
257 provenance records plus frozen runtime before/after the panel.

Both fresh subprocess fixture tests pass (isolates,self-loop,unchanged weighted
graph,import isolation,no clustering execution). Full1,283 unit tests passed
33.13s; batch syntax passed. Bounded one-CPU64GiB/two-hour batch ready for a
frozen/pushed executor before submission. No current OrthoHMM Slurm jobs remain
running; unrelated workload20915 left untouched.

Executorab364e4 pushed and frozen at publication_qfo_direct_graph_v1. Submitted
no-requeue job21326, confirmed RUNNING; output benchmarks/results/qfo_direct_graph_v1.
First worker observations pending. No scientific default or frozen runtime changed.

### Method diagram and current claim audit (2026-09-16)

Previous turn progressed via tested/frozen/submitted direct-stage diagnostic.
Objective reread;21326 confirmed RUNNING. Added method-diagram generator and
PNG/PDF/SVG manifest, grounded in the prepared factorial, full-family trace,
reconciliation trace and frozen replay source. Diagram separates initial HMM
search, cluster profiles/refinement, high-sensitivity groups, satellite candidates,
gene/species trees, reconciliation, constraints and distinct HOG/pair outputs.
It is a conceptual schematic, not execution validation. Checked frozen replay
branching: profiles start from multipass groups, not the no-profile refined
diagnostic. Constraint trace and species-tree branch shown separately; tree
bypasses and checkpoint scope disclosed. PNG visually inspected: labels/arrows
legible, no overlap/clipping. Layout test checks text bounds and pairwise overlaps.
Full1,284 unit tests passed32.96s. Manuscript links the figure.

Refreshed stale claim-checklist rows and execution gates: OrthoBench factorial,
sequence/constraint controls, corrected simulations, strata/traces and parameter
panel are complete; QfO counterparts, independent annotations, controlled scaling,
biological application and release/archive remain open. No completion claim or
expansion of scientific conclusions. Generated SVG whitespace retained unchanged
for manifest consistency.

Partial21326 progress now records minimal_imports_0 with six endpoint mismatches
before and after weights, with unchanged original constructor array; the indices
and replacements match earlier construction witnesses. frozen_imports_0 matches
at both stages. This observation does not require the optimizer or OrthoHMM imports
and precedes weight assignment, but is not an independently admitted complete
panel or a proven library/hardware cause. Remaining workers still running.

### Constructor-format diagnostic ready (2026-09-16)

Previous turn progressed with method figure and current claim checklist.
Objective reread;21326 remains RUNNING. Its first four observations now show
six pre-weight mismatches in minimal_imports_0 and frozen_imports_1, while the
other two match. Both import modes therefore have recorded failures; complete
panel admission remains pending. Constructor inputs remain intact in these rows.

Inspected installed igraph1.0.0/Python+C,64-bit IDs/abi3 extension and primary
upstream1.0.0 constructor/conversion source. NumPy enters a memoryview path;
general iterables enter a separate conversion branch. The limited-API branch
unfolds the memoryview to a list; do not assume zero-copy behavior or attribute
the failure to an unverified conditional code path. Targeted issue searches did
not establish a matching upstream explanation. Sources and caveats recorded in
QFO_CONSTRUCTOR_FORMAT_PROTOCOL_20260916.md; no installed library changed.

Extended direct-stage helper with optional Python integer-pair generator input
and fixed --compare-formats panel: three alternating minimal-import workers per
format, identical saved edge order/multiplicity/weights and vertex count. No
optimizer, partition or score. Original NumPy/default import-panel semantics
retained; running21326 uses its unchanged ab364e4 worktree. Six focused tests
pass (four fresh-worker combinations plus input-order/type and panel tests).
Full1,288 unit tests passed34.93s; batch syntax passed. One-CPU64GiB/two-hour
format batch ready to freeze/push/submit; allocations differ so clean observations
alone will not establish a fix or causal mechanism.

Executore919834 pushed/frozen at publication_qfo_constructor_formats_v1.
Submitted no-requeue21327 and confirmed RUNNING (one CPU);21326 remains
RUNNING. Outputs are separate at benchmarks/results/qfo_constructor_formats_v1.
No frozen scientific/runtime files changed and no existing job restarted.

### Direct-stage panel independently admitted (2026-09-16)

Previous turn progressed by freezing/submitting constructor-format21327. Objective
reread;21326 completed0:0 in13:27. Raw report SHA256
7b20a5e4aa440718f34d989c41ccf6e8bb427229935b6d4e1d8e051441063f8c.
Independent admission checks exact six-worker inventory, frozen ab364e4 executor,
runtime, common context, same-mode modules/libraries, import isolation,278 file
records, native observation consistency and full post-weight endpoint hashes.
All pass. Snapshot qfo_direct_graph_verified_20260916.json SHA256
08352ac4af718e5b3a9e59d1cec2d81c3a333785aeeab8f22c1ed5975cfe7ff3.

Six-edge mismatches: minimal0/frozen1/minimal2/frozen2. Clean:minimal1/frozen0.
Every stage pair has identical difference reports before/after weight assignment;
all original constructor arrays match saved inputs. Pre-weight hashes are implied
by complete bounded witnesses, not separately recorded full native hashes. Thus
weight assignment and OrthoHMM/Leiden imports are not necessary for this recorded
failure; no underlying library/hardware cause or universal historical impact is
established. No optimizer or accuracy evaluation. Manuscript, claims and dedicated
QFO_DIRECT_GRAPH_RESULTS_20260916.md updated with these limits.

18 focused admission tests passed, including incomplete/wrong-panel and corrupt
stage/hash rejection; full1,306 unit tests passed34.57s. Scientific core unchanged.
Constructor-format21327 remains RUNNING; first NumPy worker records six unchanged
pre/post-weight mismatches, other workers pending. Next admit complete format
panel and decide the next integrity-gated diagnostic from all observed outcomes.

### Constructor-format admission prepared (2026-09-16)

Previous turn progressed through independent direct-stage admission and reporting.
Objective reread;21327 remains scheduler-confirmed RUNNING. Extended the direct
audit with explicit --formats selection for job21327/executore919834/output
qfo_constructor_formats_v1. Require its exact six-worker ordered plan, minimal
imports, format agreement across plan/parent/snapshot/native result, frozen source,
runtime, preserved records and complete endpoint/hash consistency. Original21326
defaults remain intact; no partial format panel can pass.

32 focused tests passed, including cross-panel identity and constructor-format
substitution rejection; full1,320 unit tests passed35.50s. Re-admitted the actual
completed21326 report to qfo_direct_graph_admitted_recheck_v1. Its original native
report and278 checked provenance records exactly match prior admission; all six
worker summaries match after removing the newly explicit edge_format=numpy field.
Scoped diff check passes. No scientific/runtime or running-executor changes.

21327 has four worker reports so far: NumPy repeat0 has six unchanged pre/post
weight mismatches; Python-pairs0,NumPy1,Python-pairs1 match. Remaining two workers
are pending. These observations are not a completed/admitted panel or a proven
format fix. Next admit terminal results with the report hash and retain all arms.

### Format panel admitted; checked worker prepared (2026-09-16)

Previous turn progressed by tested format-panel admission support and successful
recheck of the original panel. Objective reread;21327 completed0:0 in13:53.
All three Python-pair workers match before/after weights; NumPy0 has six unchanged
mismatches,NumPy1/2 match. Raw report SHA256
d83de7fbf130642b3f4ccc5cb76618de5640a289624ab7055f38d5cd51479526.
Independent format admission passed exact plan/format/source/runtime/native
file consistency,278 file records and post-weight graph hash reconstruction.
Snapshot qfo_constructor_formats_verified_20260916.json SHA256
bfb9f49e1eb32afb31ab898d0101cc785da1b54110c70d7ea06ad0259fb87631.
Dedicated results, manuscript and claim checklist updated. Not a general fix or
causal diagnosis; no optimizer/accuracy selection or historical substitution.

Added checked_python_pair_worker.py for a subsequent bounded frozen initial-graph
replay. Requires complete admitted format panel and intact payload; records one
Python-integer-pair constructor conversion and unchanged constructor input hash.
Reuses existing full native endpoint/weight gate before optimizer and after its
return. Mismatch fails without retry. Frozen worker seed/resolution/arguments
unchanged; no scientific/default/environment edits.12 focused tests pass, including
real optimizer intact-graph execution, corrupt-graph rejection before optimization,
admission guards and restoration after constructor failure. Full1,332 unit tests
passed34.65s. Large-graph checked replays have not yet been launched; next prepare
their bounded fresh-worker orchestration and retain every partition comparison.

### Checked QfO replay executor prepared (2026-09-16)

Added run_qfo_checked_repeats.py and a no-requeue Slurm batch for exactly three
fresh sequential one-CPU initial-graph workers. Each uses the admitted Python-pair
constructor, verifies its oriented input-byte hash against saved arrays, checks
all native endpoints/weights before and after unchanged CPM optimization, checks
worker/source/runtime provenance, and validates complete unique gene coverage.
Partitions are compared without accuracy labels; disagreement is retained, not
used to retry or choose an output. Any execution/integrity failure stops the panel.
The saved graph has 976,504 vertices; this is not a full HMM replay, historical
equivalence claim, general determinism test, or controlled resource benchmark.
Unrelated GPU job 20915 was running during preparation. No existing analysis was
stopped. Full unit suite: 1,349 passed in 36.12s; focused constructor/runner tests:
29 passed, including the real optimizer's recorded arguments and full graph gate.
Next freeze this executor, submit the three-repeat diagnostic, and independently
admit its outputs before deciding how to resume the QfO ablation work.

Executor f6ad87c committed and pushed, then frozen in detached worktree
benchmarks/work/publication_qfo_checked_repeats_v1. Submitted Slurm job 21328;
squeue confirmed RUNNING with one CPU (00:07 at first check). Allocation is
64 GiB, three-hour limit, no requeue. Output:
benchmarks/results/qfo_checked_repeats_v1; scheduler log:
benchmarks/work/qfo_checked_repeats_21328.log. No result is admitted yet.
The remote still reports 21 dependency alerts (1 critical, 7 high, 11 moderate,
2 low); release/security review remains outstanding and the frozen scientific
environment was not changed.

### Independent checked-repeat admission prepared (2026-09-16)

Previous goal continuation made progress: committed/pushed the tested executor and
launched job 21328. This continuation confirmed the same job RUNNING via sacct;
no restart. At 00:04:03 the first worker's saved native-boundary record showed
all 24,148,515 canonical edge endpoints and ordered weights matching the saved
976,504-vertex graph before optimization. No post-optimizer result is claimed.

Added admit_qfo_checked_repeats.py for independent post-completion admission of
the exact job/executor/three-repeat inventory. It reconstructs graph/input hashes
using separate chunking and canonicalization, validates native before/after
observations and exact optimizer settings, verifies preserved files against the
parent report, commands, source/runtime identities and prior provenance, and
recomputes all three pairwise partition comparisons with full gene coverage.
Partition disagreement is valid evidence and is never rejected merely for being
disagreement. Admission requires a separately supplied raw-report SHA256 and
completed Slurm accounting; no incomplete run is admitted. This remains an
initial-graph diagnostic, not full HMM replay or historical equivalence.
Focused tests: 22 passed, including real igraph/Leiden observations, independent
multi-chunk hash reconstruction, invalid arrays, altered native records, exact
panel guards, and refusal to reuse output. Full suite: 1,371 passed in 35.37s.

### First matched scaling input panel prepared (2026-09-16)

Previous continuation made progress by committing/pushing independent QfO
admission2593833. This continuation confirmed job21328 remains RUNNING, including
at00:09:21, without restarting it. While that diagnostic runs, advanced the separate
practical-efficiency requirement. No completed matched scaling protocol existed.

Added MATCHED_SCALING_PROTOCOL_20260916.md and prepare_scaling_inputs.py. Basename
hash ordering with fixed salt selects nested4/8/12 complete OrthoBench proteomes;
no predictions, accuracy labels or new timing outcomes inform selection. Three
repeats each of frozen high-sensitivity, satellite_v2 and full OrthoFinder3.1.5
give27 planned runs, rotating method positions within size/repeat blocks. Planned
limits32CPUs/128GiB/24hours per run; exact command/native-runtime, simultaneous
memory/CPU accounting and throughout-run host-workload gates precede execution.
Exclusive Slurm allocation alone is not assumed to isolate non-Slurm activity.
No global cache flushing, user-process termination or timing run was performed.

Prepared benchmarks/results/publication_scaling_inputs_v1 with only symlinks to
checksum-pinned intact source FASTAs. Frozen snapshot:
publication_scaling_inputs_20260916.json SHA256
1957b050d33dd89e933ff33f96500a4fbaa6b2155d85d067d7d7d956c3d139af.
Sizes contain73,266/165,168/251,378 proteins and36,474,860/81,374,087/127,798,408
sequence characters. Independent reread checks full file inventory, input hashes,
unique coverage, strict nesting and counts. Six focused tests pass for ordering,
location independence, ambiguous input rejection, balanced27-run plan and refusal
to reuse output. Full suite: 1,377 passed in 35.00s. This is input preparation only:
taxon composition and size co-vary, no universal scaling claim, and historical
shared-node times remain descriptive rather than controlled comparison evidence.

### Checked QfO initial-graph repeats independently admitted (2026-09-17)

The interrupted prior continuation was a verified wait on live21328 and collected
resource-environment evidence: Slurm uses task/cgroup,task/affinity and cgroupv2
with ConstrainCores/ConstrainRAMSpace enabled; unrelated substantial IQ-TREE
workloads were present. No resource collector or timing run was launched. User
systemd service access exists but must not be used to move scientific jobs outside
their Slurm cgroup limits. Matched resource execution remains outstanding.

Re-read the objective and checked authoritative state: job21328 COMPLETED0:0
in38:14. Raw report SHA256
26e3bee1cc869be5df02f6160aba6fb2bb024b98a98eb348b8a5d29af12dbad5.
Independent admit_qfo_checked_repeats.py passed exact job/executor/plan/source,
runtime/command and preserved-file checks, fresh graph hash reconstruction,
before/after optimizer integrity, unique complete coverage and all three pairwise
partition comparisons.308 provenance records checked. Every output has349,898
groups and976,504 genes, byte-identical SHA256
8c162782bba955078a8f7e4e6088df5b34c9f15b055803224e513fdea8b33cdd.
Worker wall seconds1044.452333/751.549762/463.853209 are diagnostic shared-node
costs, not speedup evidence. Scheduler MaxRSS missing; no value fabricated.

Admitted output benchmarks/results/qfo_checked_repeats_admitted_v1/results.json;
snapshot qfo_checked_repeats_verified_20260917.json SHA256
78f51f5ce703caf39307c5e518ad737acdf655dbe0926a34f2c20bbdec6d03f1.
Results note, manuscript and claims updated. No accuracy evaluation, default
change, historical-output replacement or claim of general determinism. Next
extend checked construction/native-boundary observation to every payload in a
complete cached QfO replay, including multipass/profile stages, then independently
admit it before resuming QfO ablations. Scaling/application/release work remains.

### Checked full-replay payload components (2026-09-17)

Previous continuation made progress by independently admitting all three initial
graph repeats and pushing9d256a2. Re-read objective and current repository; no
OrthoHMM Slurm job remains running. Job20915 is unrelated and was left untouched.
Inspected frozen49ab replay and native isolation: exactly four clustering calls
are expected for one profile iteration, ordered initial/multipass/profile_base/
profile_expanded. All include isolates, seed4 and CPM0.1.

Added checked_replay_payload_worker.py. It requires the admitted repeat snapshot,
checks its308 prior provenance records, validates a separately hashed per-stage
payload manifest, identical QfO gene order, metadata, array types/shapes/endpoints
and finite weights, and refuses already observed payloads. It uses the existing
Python-pair adapter plus frozen worker/native before-after graph observer with
one-CPU affinity, recording helper/source/admission provenance before execution.
It does not require later-stage edges to equal the initial graph: those edges
must match their own freshly preserved payloads.

Added checked_replay_interceptor.py to intercept only the exact frozen isolated
worker invocation. Original graph generation and serialization remain unchanged.
Copies all five payload files before the temporary directory disappears, checks
original/copied bytes, records commands/manifests/logs/failures, requires a caller-
provided result validator before continuing, and preserves each returned partition.
It rejects a fifth call or any continuation after a failed call; unrelated
subprocess calls pass through unchanged. No existing frozen files were edited.

Ten focused tests pass for stage/input/settings guards, duplicate-observation
rejection, preservation, failure/no-retry behavior, ordinary subprocess forwarding,
and exact four-stage inventory. Synthetic interceptor tests do not establish
native full-replay success. Full suite: 1,387 passed in 34.93s. The parent launcher,
per-stage post-worker validator and independent completed-replay admission still
need integration before submitting a full cached replay. No new inference or
accuracy evaluation was launched in this continuation.

### Complete checked QfO replay launcher prepared (2026-09-17)

Previous continuation made progress by committing/pushing payload components
80a0166. Re-read objective and inspected current sources. Added
validate_checked_replay_payload.py to verify actual worker modules/environment/
metadata, admission/helper manifests, original oriented constructor bytes and
complete native endpoint/weight hashes before/after unchanged optimization.
It validates unique complete partition coverage before the replay consumes it.
Initial-worker preflight now also requires exact previously admitted graph bytes;
later graphs remain bound to their own generated payloads.

Added run_qfo_checked_full_replay.py. Parent verifies prior308 provenance records,
frozen7f/49ab native runtime and input audit before running; a fresh worker imports
the frozen replay first, then intercepts only its isolated clustering calls.
Exactly four successful checked calls, four full-coverage final stage partitions
and nonzero constructed HMM profiles are required. Input/runtime/source checks
repeat after execution. Initial partition comparison is recorded without rejecting
disagreement or selecting by accuracy. Original payloads, failures and returned
partitions remain preserved; no retries or scientific/default changes.

Prepared no-requeue batch:32CPUs for profile work, one CPU per native graph
worker,192GiB,24-hour limit. Disk inspection shows approximately12TiB available.
This is a shared-machine cached replay with observation overhead, not controlled
end-to-end efficiency evidence. Nine new validator/parent tests plus ten existing
interceptor tests pass; validator fixtures use actual igraph/Leiden observations
and test changed graph hashes, constructor hashes, metadata, scientific module
records, provenance and incomplete/duplicate output membership. Full suite:
1,396 passed in35.04s. Freeze/push executor before submitting; independent
completed-replay admission remains required before QfO ablation inputs are used.

Executor96333fd committed/pushed and frozen in detached worktree
benchmarks/work/publication_qfo_checked_full_replay_v1. Submitted job21329;
squeue confirmed RUNNING with32CPUs at00:07. Output directory:
benchmarks/results/qfo_checked_full_replay_v1; scheduler log:
benchmarks/work/qfo_checked_full_replay_21329.log. No replay output admitted yet.
The existing21 remote dependency alerts remain outstanding for release review;
this submission did not change the frozen dependency environment.

### Full checked-replay independent admission prepared (2026-09-17)

Previous continuation made progress: connected/tested/froze executor96333fd,
submitted21329 and pushed its ledger recordc69b81e. Re-read objective and polled
the same handle: RUNNING at00:04:17. The preserved initial-stage native record
shows the admitted complete endpoint hash before optimization; no completed
stage or full replay result is inferred from that observation. No restart.

Added admit_qfo_checked_full_replay.py. Requires exact completed job21329 and
executor96333fd, the complete four-call/four-stage inventory, unchanged full
replay parameters, nonzero built profiles, source and command inventories,
preserved-stage execution/provenance, all-stage gene order, original/copied
payload hash agreement, independent graph reconstruction and native before/after
integrity. Recounts every retained partition and final-stage coverage, verifies
that multipass/profiles outputs came from their corresponding checked workers,
and crosschecks recorded initial/multipass edge counts. Re-runs the frozen input
auditor and native-runtime verification before producing admitted results.

The historical final partition is compared against profiles_refined and reported
without treating disagreement as a failed experiment or selecting by accuracy.
Temporary original payload paths are retained as provenance but not incorrectly
required to still exist. Every permanent retained payload and reported provenance
record is rechecked. Exact raw-report SHA256 must be supplied after completion;
an existing admission directory is never reused. Fifteen inventory tests pass,
covering incomplete/failed/wrong-job/wrong-order/changed-parameter/scored runs,
absent HMM profiles and existing-output rejection. Full suite:1,411 passed in35.75s.
No output has yet been admitted, and no benchmark score/default has changed.

### Read-only Slurm resource-accounting feasibility (2026-09-17)

Previous continuation progressed independent full-replay admission7418845.
Re-read objective and confirmed21329 still RUNNING at00:06:11, same initial
clustering stage. Advanced the separate practical-efficiency requirement while
it runs. scontrol listpids identified anchor3587494 and two descendants inside
/system.slice/slurmstepd.scope/job_21329/step_batch/user/task_0. No process moved,
limit changed, cgroup peak reset, or unrelated workload stopped.

Added slurm_resource_snapshot.py: validates exact requested job/step scope and
unified cgroup, retains raw/keyed cpu.stat, memory.current/peak/stat/events,
effective cpuset and memory limits through job ancestry, and samples descendant
RSS with process identity/membership checks and explicit race/error records.
CPU is cumulative microseconds; cgroup peak is since creation or last reset, not
an inference-only measurement. Summed RSS is non-atomic and may double-count
shared pages; cache/kernel cgroup charges are not RSS. Host isolation is not
established by this collector. Updated matched-scaling protocol accordingly.

Live snapshot at2026-09-17T14:39:33.689029+00:00 retained three processes, no
sampling errors,32 effective cgroup CPUs and192GiB inherited memory limits.
Reported cgroup peak6,886,125,568 bytes; summed sampled RSS6,456,483,840 bytes;
cumulative CPU575,346,781 microseconds. These are partial-run diagnostic values,
not final costs or speedup evidence. Snapshot:
benchmark_tools/results/qfo_checked_full_resource_snapshot_20260917.json.
Original observation preserved separately; v2 corrects peak terminology to
allow creation-or-reset rather than asserting an unverified reset history.

Nineteen focused tests pass for scope/path guards, malformed counters/cpusets,
raw preservation and distinct memory quantities. Full1,430 tests passed36.64s;
nineteen focused tests rerun after the peak-field terminology correction.
Time-series monitoring, controlled workload gates, stage timing, failure/timeout
accounting and actual27-run scaling execution remain outstanding. No new
scientific inference was launched and no accuracy/default changed.

### Bounded resource time series and replay stage progress (2026-09-17)

Previous continuation progressed read-only cgroup/RSS accounting860cb0e.
Re-read objective and confirmed21329 RUNNING at00:10:55. Later inspection found
cluster_0_initial/execution.json statuschecked,349,898 groups/976,504 genes and
partitionSHA8c162782bba955078a8f7e4e6088df5b34c9f15b055803224e513fdea8b33cdd,
matching admitted initial repeats. cluster_1_multipass is now present and running.
This is stage-level progress, not completed replay or independent full admission.

Added monitor_slurm_resources.py for bounded read-only JSONL telemetry. It
checks stable anchor PID creation time/job/scope and nondecreasing CPU counters,
retains each observation before proceeding, summarizes CPU deltas over actual
sample span, keeps reported cgroup lifetime/reset peak separate from sampled RSS
maximum, and retains failures without asserting that the scientific job stopped.
Host load and per-CPU counters are recorded but do not attribute unrelated load
or certify exclusive execution. Existing output is refused; nothing is restarted.

Live smoke used confirmed anchor3587494/job21329, five samples at2-second target
interval. All five retained, span8.003259279s, CPUdelta8,000,914microseconds,
sampledRSSmaximum6,861,684,736bytes, reportedcgrouppeak6,886,125,568bytes, no sampled
process errors. These are partial-run diagnostic values, not final inference costs.
Unit tests and unrelated jobs were running; no matched-speed claim is made.
Reports qfo_resource_series_smoke_20260917.json and
qfo_resource_series_samples_20260917.jsonl copied to benchmark_tools/results.
Independent reread verified sample hash/count, CPU delta, RSS maximum and stable
anchor scope. Raw outputs preserved in benchmarks/results/qfo_resource_series_smoke_20260917.

Fourteen focused tests pass for exact deltas, incomparable identity/counter/time
rejection, retained partial failures, successful collection, invalid plans and
existing-output refusal. Full1,444 tests passed37.13s. Updated matched-scaling
protocol; stage-boundary timing, whole-run monitoring, unrelated-workload gates,
timeout/failure accounting and actual27-run scaling execution remain outstanding.
No new inference, benchmark score or default change in this continuation.

### Measured command lifecycle wrapper (2026-09-17)

Previous continuation progressed telemetry78cc9e1. Re-read objective;21329 remains
RUNNING at00:15:24. Added measure_slurm_command.py to measure one fresh command
inside its existing dedicated Slurm task subtree. Requires only the wrapper in
that subtree before launch, exact effective CPU count and inherited task/job
memory cap. Retains command/log/exit status and samples across launch-to-exit,
distinguishes zero/nonzero exit, timeout, spawn failure and measurement failure.
Sampling failure does not terminate a still-running scientific command: native
execution continues to exit or the declared timeout, with invalid measurement
flagged. Only its newly created process group is targeted for timeout cleanup;
TERM-resistant children receive KILL after a grace period. No Slurm job/cgroup
is moved or reconfigured, and no unrelated process group is signaled.

Twelve initial tests with real subprocesses passed; full1,456 tests passed36.91s.
Then added a real parent/child TERM-resistance test: all13 focused tests pass.
The wrapper includes its own cgroup accounting overhead and cannot establish
host quietness, complete biological output or inference-only cgroup memory peak.
Prepared a two-CPU/one-GiB no-requeue Slurm smoke test with small parent/child
memory allocations and a brief child calculation/sleep. Bash syntax passes.
Freeze/push before execution. This smoke is instrumentation validation, not one
of the planned27 scientific scaling runs, and does not count as timing evidence
for any orthology tool. Controlled workload gates and exact scientific command/
output validation remain outstanding.

Executor8c54f87 committed/pushed and frozen at
benchmarks/work/publication_resource_command_v1. Submitted instrumentation
smoke21330; authoritative sacct confirms COMPLETED0:0 in00:03. Wrapper status
command_exited_zero, exit0, command wall2.424314534s;13 samples spanning2.427225855s.
Observed process counts were1 before launch,3 during execution,1 after exit;
two effective CPUs throughout and inherited1GiB cap passed preflight. No process
sampling errors. Cgroup CPUdelta450,593microseconds; reported memory peak69,480,448
bytes and sampled summedRSSmaximum84,934,656bytes retain their different accounting
definitions (shared pages can make summed RSS exceed cgroup charge).

Independent reread checked retained sample hash, scope CPU count, before/during/
after process inventory and cumulative CPU difference against the report.
Results and raw samples copied to resource_command_smoke_20260917.json and
resource_command_smoke_samples_20260917.jsonl under benchmark_tools/results;
originals remain at benchmarks/results/resource_command_smoke_v1. No quiet-host
claim: controlled_workload_verified staysfalse. Updated matched-scaling protocol.
QfO21329 remains on multipass clustering; its final replay admission is pending.

### Simulation generating-tree controls prepared (2026-09-17)

Previous continuation progressed measured-command wrapper8c54f87 and successful
Slurm smoke4c22d15. Re-read objective and confirmed21329 still RUNNING at00:22:53
and00:28:56. Advanced the separate simulation-truth tree-error requirement while
the replay runs, rather than selecting cases by already inspected accuracy.

Added SIMULATION_TREE_CONTROL_PROTOCOL_20260917.md and
prepare_simulation_tree_controls.py. All70 fixed variable-length condition/seed
cells undergo fresh generation/history/input checks using manifest806aa1e5...
The parent generating ExtantTree is pruned to each dataset's retained taxa,
checking unique names/binary topology/nonnegative lengths and preserved pairwise
patristic distances. Existing tested rooted-NNI helper supplies a deterministic
generating/one-NNI/two-NNI panel at rooted clade distances0/2/4. Perturbed lengths
travel with subtrees; these are oracle/stress diagnostics, not empirical tree
uncertainty or achievable end-to-end information.

All70 cells prepared successfully with210 trees and420 unique planned runs for
satellite_v2/fullOrthoFinder3.1.5, without filtering originally failed methods.
Protocol requires exact CLI/root/name/output semantics, unchanged supplied-tree
mode controls and upstream identity before inference/scoring attribution. Paired
seed-level mean differences use20,000 resamples, seed20260918 and126 exploratory
F1/P/R endpoints, retaining failures and multiplicity slots. No new defaults or
held-out-validation claim; original variable-length outcomes are development-exposed.

Prepared benchmarks/results/simulation_tree_controls_prepared_v1; snapshot
simulation_tree_controls_prepared_20260917.json SHA256
0caffab16f73019fabafe1fcfaac8c2de90960c5f1317803bf83308a6d2bf914.
Independent reread checked all210 tree hashes, unique exact taxa and clade
distances against original generating-tree clades projected onto retained taxa,
plus420 distinct condition/seed/tree/method combinations. Six focused tests pass;
full1,463 unit tests passed37.60s. Claims updated to distinguish preparation from
execution. No supplied-tree method run or score has yet been produced.

### Simulation native tree-format compatibility (2026-09-17)

Re-read the objective. The immediately preceding prompt-only response was not
scientific progress; revalidated repository/scheduler state before continuing.
Confirmed73514f0 local and remote. QfO21329 remains RUNNING at00:35:45 and00:42:02;
no restart or duplicate replay was submitted.

Audited installed OrthoFinder3.1.5 supplied-tree validation/rooting and frozen
OrthoHMM supplied-tree parsing. Native parser testing exposed a real blocker:
the prepared `[&R]` prefix makes OrthoFinder see a NoName tree and reject all
expected taxa on the first input. OrthoHMM accepted all210 original files.
Preserved the frozen original manifest/files. Added a structured plain-Newick
adapter requiring exact descendant-clade/branch-length identity after reread.
Generated210 derivative trees; both native parsers then accepted all210.
These are input-compatibility checks, not reconciled method results.

Derivative snapshot simulation_portable_trees_prepared_20260917.json SHA256
b9ed4fb8dc27da28dd56c674d1ece2edbb3a04697dec14ea3d4538bf6d9dbc0b.
Original generated output benchmarks/results/simulation_portable_trees_v1.
Preparation status remains pending native checks because direct terminal checks
are documented evidence rather than a standalone machine-readable admission.
Execution preflight must repeat them and verify all tree/source hashes.

Added fresh_supplied_method command builder preserving frozen non-tree arguments
and refusing existing outputs, overlapping paths and restart/pre-supplied flags.
Dry-constructed420 unique fresh-run destinations from pinned source manifests;
no inference launched. Actual commands must reference the derivative trees.
Local source shows supplied OrthoFinder trees bypass STRIDE/multiple-root
handling, reinforcing the unchanged-tree mode-control requirement. Cached -ft
reuse and fresh-run upstream identity are not yet validated. Do not attribute
prediction changes to topology before those gates pass.

Eighteen new focused tests passed; full unit suite1,481 passed36.52s. Updated
simulation protocol with parser incompatibility, derivative provenance and
execution semantics. Next: freeze and execute unchanged-tree mode controls with
native source/input/output admission before the420 main method/tree runs.

### Unchanged-tree execution pilot prepared (2026-09-17)

Previous continuation progressed committed26b97ca native parser compatibility.
Re-read objective and confirmed21329 RUNNING at00:44:01, then observed its
authoritative terminal FAILED1:0 at00:46:36. It completed checked initial and
multipass clustering; profile_base native execution returned but the validator
rejected OMP_NUM_THREADS=32 instead of1. Frozen profile_expansion.py:702 sets
this variable in the parent. Failure and outputs are preserved without retry;
the child environment needs explicit isolation before any corrected replay.

Added run_simulation_tree_mode_control.py and a four-CPU16GiB one-hour batch
for baseline_20261101, the first baseline seed, as an execution pilot, not a
selected accuracy result. Both methods rerun fresh using their own previously
inferred tree. No cache reuse, source changes or accuracy scoring. Native
admission, pair membership, rooted topology and retained artifact hashes are
reported separately; postprocessed OF artifacts alone do not prove every
upstream state. Independent admission and all-dataset mode controls remain.

Preflight freshly re-admitted both original methods, verified806inputgenes,
and passed both native tree parsers. Retained inventories contain27OH and202OF
artifacts. New runner preserves failures and refuses existing destinations.
Nine focused tests pass; full1,490unit tests passed36.49s. Executor will be
committed/frozen before scheduler submission; no pilot results claimed yet.

### Mode pilot complete; QfO child-environment correction (2026-09-17)

Committed/pushed0e797ba and froze publication_simulation_mode_control_v1.
Initial submission21331 failed its exact-commit shell guard after1s because the
submitted commit argument was mistyped; no inference output was created.
Confirmed terminal failure, corrected the argument, and submitted21332.
Job21332 COMPLETED0:0 in35s. Both fresh native runs admitted by the existing
validators and retained identical rooted topologies and native ortholog pairs:
OrthoHMM2,881pairs; OrthoFinder2,884pairs. No accuracy calculation performed.
OrthoHMM's27retained artifacts matched byte-for-byte. OF retained201versus202:
the species-tree alignment is absent in supplied mode and the two MCL files
differ in command-comment paths (inspected diff); independent semantic reread
remains required. Pilot is not all-dataset equivalence or publication admission.
Snapshot simulation_mode_control_baseline_seed1_20260917.json retains pending
independent-admission status and all comparisons, including mismatches.

Preserved QfO21329 failed parent/worker snapshots and documented the exact
environment mismatch in QFO_CHECKED_REPLAY_ENVIRONMENT_FAILURE_20260917.md.
All three observed native boundaries report intact graphs; profile_base is
still rejected, not silently admitted. Explicit child-only OMP/OPENBLAS/MKL=1
overrides now prevent profile-stage parent settings from leaking into clustering.
All other environment values are inherited; parent/profile settings unchanged.
New four-call regression confirms this isolation and records inherited values.
Twenty focused clustering-wrapper/validator tests pass; full1,491unit tests
passed35.90s. Prepared fresh v2 batch; freeze/submit after committing. Existing
independent full-replay admission remains pinned to v1 and needs an explicit
v2 provenance update before any completed retry can be used.

### Corrected replay submitted and admission pinned (2026-09-17)

Committed/pushedc93eb2c and froze publication_qfo_checked_full_replay_v2 at
c93eb2c8cf5671b9e99b3e534f68d11fac6282d7. Submitted21333 using the commit read
directly from git; confirmedRUNNING00:00:12. Original input audit passed
976,504genes78FASTAs, with the same7documented historical source differences.
Output benchmarks/results/qfo_checked_full_replay_v2; scheduler log
benchmarks/work/qfo_checked_full_replay_v2_21333.log. Do not restart while live.

Extended independent full-replay admission with explicit v1/v2 choice and fixed
job/revision/output mappings. V2 requires recorded inheritedOMP1/1/32/32 and
childOMP/OPENBLAS/MKL1 for all four stages, alongside unchanged native-worker
environment/graph/source/coverage gates. Six additional tests cover version/job
and inherited/override mismatches; all21focused admission tests pass. No v2
result is admitted yet.

Independent OrthoFinder checkpoint reread for pilot21332 recovered98groups and
806genes on each side with identical canonical memberships, confirming the
observed MCL byte differences do not change this pilot partition. Added bounded
pilot summary with snapshotSHA446a3c19cebf87fe562cc9e88cf183986e4e644b71116eaa5ba2839c9083f561.
Complete independent pilot admission and the remaining dataset controls are
still outstanding; no general equivalence or accuracy advantage claimed.

### Pilot independently admitted; all-dataset mode controls prepared (2026-09-17)

Previous continuation made progress via0e797ba/c93eb2c/3c853f4: pilotexecution,
QfOdiagnosis/correction and live retry. Re-read fullobjective and confirmed
21333RUNNING00:03:11; no duplicate replay. Added fixed pilotadmission script
rechecking21332scheduler,0e797baexecutor/sources, input/output inventories,
native baselineadmissions, independently reconstructed commands and all pairs,
rootedtopologies and retainedartifacts. Both methods passed. Exact exceptions
allow the absent OFspecies-tree alignment and one command-comment line per MCL
file, not any matrix change. Nine focused tests cover these gates.

Independent snapshot simulation_mode_pilot_verified_20260917.json SHA256
4a4630c8c036a4ae85045045ebf848ab9eefa730de9dcc8295570a50deeddd0c.
No accuracy scoring; one-cell equivalence does not imply general equivalence.

Prepared simulation_mode_panel_prepared_20260917.json SHA256
03f0d115d0f2a1f5354ca6a8a369e419dbd9e10945e450dfd2534d61212bb747:
all70datasets/140methodslots,130newruns,2pilot reused,8unavailable original
baselines (3OH/5OF). The runner now accepts a canonical subset of methods so one
original failure cannot suppress the other tool. Panel selection uses pinned
completion states only and preserves every failure reason; all8remain in future
420oracle/NNI run inventory. Added per-task evidence/failure retention and
4CPU16GiB array with max2concurrenttasks. Thirteen panel/subset tests plus9pilot
and9existingrunner tests pass(31total). Freeze/submit after full tests/commit;
no all-dataset mode equivalence or completed oracle experiment claimed.

Fullunit suite1,519passed35.95s before executor freeze.
