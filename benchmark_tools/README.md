# Benchmark Workflows

## HMM Evidence Campaign

The [campaign protocol](hmm_evidence_campaign_20261018.json) and
[execution ledger](results/hmm_evidence_campaign_status_20261018.json) describe
the bounded investigation of weighted profiles, alignment-aware recruitment,
local joint-tree rescue and selective HMMER validation. This is development
work, not a new validated default or a completed publication release.

Weighting and the standalone alignment guard have production opt-in flags
described in the [main README](../README.md). Local rescue and selective
HMMER are experimental package/workflow integrations; promising results
require complete-dataset follow-up and a tested production opt-in interface
before recommendation. HMM similarity alone does not establish orthology.

The [completed weighting/local comparison](results/hmm_completed_method_comparison_20261019.tsv)
and [readback summary](results/completed_weighting_local_readback_20261019.json)
retain all eight endpoints and the contextual tools. Both weighting variants
pass the original development screen but have mixed primary-benchmark effects;
neither is a default recommendation. Local rescue leaves all 70 reference-family
scores unchanged and stops at its original gate. The failed original HMMER panel
is explicitly unadmitted, not a negative accuracy result. This interim table does
not replace the required controlled recovery or final campaign report.

The first selective HMMER biological pilot failed a source-residue integrity
check, rather than producing an admissible negative accuracy result. The
[projection recovery protocol](hmm_projection_integrity_recovery_20261018.json)
records the diagnosed center-star insertion-projection fault and a prospective
matched recovery panel. Package-level `projection_policy="residue_preserving"`
is an opt-in repair in `center_star_msa` and `build_msa_profile`; the default
remains `legacy`, and frozen scientific jobs are unchanged. HMMER reconstruction
uses the profile's recorded policy without relaxing source/profile checks.
The repair is also available as experimental API `profile_projection_policy`
and CLI/harness/replay `--profile-projection-policy residue_preserving`.
Production requires fresh built-in high-sensitivity inference; replay requires
the original sequences. Worker and jackknife profiles retain the policy.
Untested combinations with weighting or alignment guards are rejected. Metrics
identify the selected policy, and production records the repaired-profile count.
The five-arm recovery runner is tested, including a complete small native
fixture. Its [independent reader contract](hmm_projection_recovery_readback_20261019.json)
requires separate repair/decision/backend contrasts and exact repaired-control
identity. `report_projection_recovery.py` binds an explicit full immutable launch
revision, the original parent job and complete five-arm producer, and repeats
native/score/ownership readback before admission. Freeze committed/pushed runner
and reader sources before biological submission; pass the new immutable launch
revision to the separately frozen reader. The original failed HMMER/collator
bindings are never redirected. This is not a validated accuracy improvement or
a new production recommendation.

The [Stockholm annotation recovery](hmm_stockholm_annotation_recovery_20261019.json)
addresses a separate post-build validation error: HMMER's annotated re-save can
represent fragment-end missing data with `~` and remove empty insertion columns.
Validation now preserves exact source ownership, every residue including `X`,
and model match-column count while interpreting documented nonresidue symbols.
The independent reader repeats this check from raw alignments. Construction,
scoring, fragment settings and original advancement gates remain unchanged.
Use a fresh complete five-arm pilot and explicit reader
`--launch-entry stockholm_annotation_recovery_pilot`; do not retarget or admit
the failed original recovery or its dependent jobs.

The [explicit Stockholm transfer binding](hmm_stockholm_projection_followup_20261019.json)
uses `plan_projection_followup.py --campaign stockholm_annotation` with immutable
science24519/reader24520 launch snapshots. The default planner still targets the
original recovery. Both enforce the same full-panel admission, selection gates,
inactivity rules and minimal matched-control closure. The new option prepares
inputs only after successful independent readback; it does not redirect the old
failed preparation or authorize the old transfer driver. The standalone current
driver accepts `prepare --campaign stockholm_annotation` and reauthenticates the
new frozen planner/source and preparation job24521 in an isolated process. It
preserves all original scientific settings, control dependencies and exact
measurement-assessment reuse. Its default remains the original recovery.
Separately frozen, tested and committed/pushed execution and full-vector reader
bindings remain required before new QfO/Three Kingdoms inference.
The [new full-vector reader contract](hmm_stockholm_projection_followup_readback_20261019.json)
uses `report_projection_followup.py --campaign stockholm_annotation
--launch-revision <full immutable commit>`. It requires the explicitly named new
execution entry with both submitted dataset jobs, exact new preparation and
separate reader bindings. Default reader behavior remains the original campaign;
old/partial/cross-campaign evidence cannot enter the new report. Native/score,
functional-population, matched-control and sequential-accounting checks are
unchanged, and this report alone cannot close the campaign or promote defaults.

| Workflow | Execution | Separate result validation |
| --- | --- | --- |
| Weighted-profile and alignment-guard OrthoBench pilots | `run_profile_weighting_pilot.py` | `report_profile_evidence_pilot.py` |
| Unchanged weighting transfer to QfO and Three Kingdoms | `run_profile_weighting_followup.py` | `report_weighting_followup.py` |
| Local nomination, cached joint-tree readback and final native pilot | `run_local_rescue_probe.py`, then `run_local_joint_pilot.py` | `report_local_joint_pilot.py` |
| Matched off/builtin/HMMER-measurement/HMMER-reselect pilot | `run_profile_weighting_pilot.py prepare --experiment selective_hmmer`, then run its plan | `report_profile_evidence_pilot.py --experiment selective_hmmer` |
| Matched five-arm projection/HMMER recovery | `run_profile_weighting_pilot.py prepare --experiment projection_recovery`, then run its frozen plan | `report_projection_recovery.py --launch-revision <full immutable commit>` |
| Combined admitted results and explicit follow-up decisions | No new inference | `report_hmm_evidence_campaign.py` |
| Separately bound recovery campaign comparison | No new inference | `report_hmm_recovery_campaign.py` |
| Gated projection transfer input preparation | `plan_projection_followup.py` after independent pilot admission | Inputs only; separately frozen execution/reader and original gates required |
| Fixed projection transfer execution | Standalone `run_projection_followup.py prepare`, then its frozen Slurm driver | `report_projection_followup.py --launch-revision <full immutable commit>` |

The [recovery report protocol](hmm_recovery_campaign_report_20261019.json)
retains the original successful weighting/local reader bindings and separately
binds the complete five-arm recovery reader. The original failed HMMER panel
remains an explicitly unadmitted all-NA row. Repair-only gains, matched decision
effects and HMMER-minus-builtin effects remain distinct. Functional-population
comparisons qualify FAS eligible pools and carry the previously audited numeric
sample counts; these are conditional unseeded means, not orthology F1 or paired
independent observations. This report does not close the campaign or promote
defaults while follow-up/confirmation requirements remain.

The [projection transfer contract](hmm_projection_followup_20261019.json)
requires complete original pilot/reader success and unchanged gates before
preparation. It retains every passing candidate and the original inactivity
fallback, with matched construction/decision/measurement controls explicitly
separated from advancement. A failed-gate arm can be required as a control but
cannot thereby become an advancing candidate. If nothing qualifies, preparation
records no follow-up rather than searching additional datasets for a pass.

The projection execution driver reauthenticates the original successful
preparation job and independently admitted pilot in a separate frozen-planner
process. It reconstructs the complete original input/source/cache inventory,
preserves the original scientific namespace, and enforces matched-control
dependencies even when a control is also an advancing candidate. Measurement-only
arms must reproduce repaired graph/candidate/native/tree products and reuse the
exact repaired assessment with fresh score readback; QfO FAS is not resampled for
this control. Preparation and execution preserve previous attempts and use short,
distinct fresh scoring identifiers. This driver does not replace an independent
post-run reader, authorize a new biological submission, or promote defaults.

The [independent projection transfer reader](hmm_projection_followup_readback_20261019.json)
binds both actual dataset jobs, exact execution plans and separately frozen reader
identity through an immutable launch snapshot. It verifies complete ownership,
commands, raw resources, profile/backend evidence, native products, control
equivalence and source/runtime fingerprints before admitting eight-endpoint
vectors. Scheduler accounting must show nonoverlapping QfO-before-Three-Kingdoms
inference; this is not a quiet-host timing requirement. Newly audited raw
functional populations are content-pinned against the admitted aggregates,
including numeric FAS sample counts distinct from the eligible pool. Stopped,
applicability and measurement roles remain explicit; repair-only gains cannot
establish backend benefit. The reader is not independent biological confirmation
and does not itself authorize default promotion or goal completion.

Use the ledger's immutable execution plans and frozen source revisions rather
than rerunning these scripts with current modules or reusing an output folder.
The weighting follow-up driver must run as a standalone file, not with `-m`;
it binds the original pilot's scientific modules before importing benchmark
packages. Preparation, inference and result validation have distinct commands;
consult each script's `--help` and the retained exact submission records.
Weighting and alignment pilots likewise select `prepare --experiment weighting`
or `prepare --experiment alignment` before running their frozen plans.

The post-run readers require complete successful source jobs, raw content and
score checks, and the original controls/gates. Existing cached files are not
proof of completed inference. The combined report retains all eight endpoint
columns, explicit `NA` for unrun endpoints, precision/recall/coverage and
resource details in its JSON, and conditional EC/GO/FAS assessed populations.
Functional counts preserve the native `NR_ORTHOLOGS` coordinate. For FAS this
is an eligible-pair pool before subsampling, not the number of numeric scores
in the saved mean; distinguish those counts when interpreting coverage or
uncertainty. The [real-control helper check](results/weighting_qfo_control_helper_integration_20261018.json)
records both for the canonical QfO control without creating a new FAS sample.
Measurement-only HMMER must reproduce off and cannot advance. Failed fixed
branches remain stopped without retuning; an unmeasured pool requires its
prespecified applicability fallback, not a negative accuracy claim.

OrthoBench and corrected QfO are primary; TreeFam is supporting and Three
Kingdoms supplementary. Existing panels are development-exposed. Functional
means are not F1, FAS is unseeded, and shared-Threadripper cached replay timing
is approximate, not isolated fresh-pipeline timing or aggregate process-tree
memory. Passing a point screen does not establish significance or generalize
to arbitrary datasets. The collator does not promote defaults or close the
campaign; its `required_next_actions` remain work to complete.

## Publication Study

Start with the [publication reproduction guide](PUBLICATION_REPRODUCTION.md)
for the frozen method, corrected QfO release, admitted results, statistical
reproduction and outstanding release gates. The manuscript is a working draft,
not a completed publication release.

The [private Threadripper executor route](THREADRIPPER_PRIVATE_EXECUTOR.md)
now explicitly binds the retained private plan, lookup and controller without
falling back to shared Python. All 229 focused checks pass; no readiness permit
or timing run was created, and controlled timing remains deferred.

The [latest nine-page main review](results/PUBLICATION_MAIN_SIMULATION_REVIEW_20261002.md)
now includes the completed evolutionary simulations and supplied/perturbed-tree
results, preserving negative findings and failures. All pages were inspected;
earlier renders remain historical and no scientific scores or defaults changed.

The [standalone simulation arithmetic replay](SIMULATION_ARITHMETIC_REPLAY.md)
checks both frozen length panels from their existing reports after copying
outside the checkout. It preserves fixed-length unavailable comparisons and
variable-length conditional contrasts; no native inference or scoring is repeated.

The [standalone YGOB arithmetic replay](results/YGOB_ARITHMETIC_REPLAY_RESULT_20261002.md)
now checks all four methods' point metrics and six paired interval endpoints
from identifier-free retained counts after fresh local archive extraction.
Scores stay frozen; this is not native re-admission or independent-family
validation. Preceding public CI at `a485ba00` has all eight jobs successful,
with its four raw-source exclusions and remaining execution-time skips explicit.

The [public and raw-data test profiles](results/PUBLIC_TEST_PROFILE_20261002.md)
separate four explicitly declared SwissTrees raw-source cases from public CI,
without changing the default full suite or exporter assertions. The collection
gate rejects incomplete inventories; all 72 focused cases pass with real raw
inputs. A public-profile pass is not full raw/native reproduction or release readiness.

The [current accuracy overview](results/CURRENT_ACCURACY_OVERVIEW_20261002.md)
now matches the unchanged corrected score table and contemporary matched-input
SonicParanoid Three Kingdoms result. It verifies all 72 table cells and 48
corrected QfO coordinates without native rescoring; input-consumption gaps and
metric/output distinctions remain explicit. Older figures stay historical.

The [prespecified synthetic null-score diagnostic](results/FROZEN_NULL_SCORE_RESULT_20261002.md)
retains 90,000 independent random pairs and all 90 tail endpoints. Strong
glutamine bias frequently passes the approximate significance filter; ordinary
background samples cannot establish rare-tail calibration. This is not real-data
orthology accuracy, a default change or controlled timing. All scores remain frozen.

The [latest eight-page main review](results/PUBLICATION_MAIN_NULL_REVIEW_20261002.md)
adds the prespecified null-score protocol, findings and significance limitation.
All pages and the new nine-panel figure were inspected; scientific scores and
defaults stay frozen. The [earlier seven-page review](results/PUBLICATION_MAIN_REVIEW_20261002.md)
preserves its descriptive-table/private-input restoration snapshot. Neither
is submission-ready or a new complete study archive.

The [earlier seven-page manuscript review archive](results/PUBLICATION_REVIEW_COMPONENT_RESULT_20261002.md)
now verifies after fresh extraction outside the checkout: 55 payload files and
34 direct local targets. It preserves the committed review and its render-time
ledger, not current CI status. Transitive evidence, complete native workflows,
raw-data rights and public release remain outside this review-only component.

The [native CI utility follow-up](results/CI_NATIVE_UTILITIES_20261002.md)
records bounded macOS confirmation of fresh CPU kernels and cgroup/tmpfs test
fixtures. It adds explicit GNU Time/Pandoc test setup with 100 passing local
checks, not production timing, full remote CI success or release clearance.
Frozen scientific commands and resource admission remain unchanged.

The [Three Kingdoms assessment fixtures](results/CI_THREE_KINGDOMS_REFERENCE_FIXTURES_20261002.md)
exercise real conversion/scoring with separate synthetic and pinned BUSCO
reference branches. Missing biological data skips only its labeled branch;
synthetic tests do not replace production accuracy or clear data rights.
It also removes fixed true/false paths from three QfO smoke fixtures. All 158
combined local cases pass. The subsequent actual macOS fast log confirms 152
passes and six retained-data skips within that scope, not full-suite success.

The [launcher environment follow-up](results/CI_LAUNCHER_ENVIRONMENT_20261002.md)
contains test-only environment leakage and relocates a temporary replay fixture.
All 94 ordered local cases pass without environment deltas. New macOS confirmation
of this fix remains open; historical commands and scientific results stay unchanged.

The [additional launcher follow-up](results/CI_ADDITIONAL_LAUNCHER_ENVIRONMENT_20261002.md)
contains six more passing-but-leaking synthetic launch cases using an explicit,
non-autouse fixture. All 183 ordered local cases pass without environment deltas.
The earlier replay fix is macOS-confirmed; complete CI/release remains unfinished.

The [retained-record follow-up](results/CI_RETAINED_RECORD_BINDINGS_20261002.md)
uses explicit checkout test bindings without repinning scores/manifests, and
strengthens the standalone OrthoBench source-identity guard. All 85 local cases
also pass in a copied tree with original-checkout reads blocked. This does not
provide a complete portable native workflow or close remaining release gates.

The [owned-group timeout follow-up](results/CI_OWNED_GROUP_TIMEOUT_20261002.md)
shares bounded subprocess cleanup between two wrappers and requires proof that
a group is absent before accepting a permission-error race. All 55 local cases
pass; new macOS confirmation is open. Prospective inventories must pin the new
helper, while historical receipts and scientific results remain unchanged.

The [exit-waitability follow-up](results/OWNED_GROUP_WAITABILITY_20261002.md)
addresses a later EPERM before the owned leader becomes waitable. A bounded
polling window still requires fresh absent-group proof; live or denied groups
remain failures. All 135 focused local cases pass. New macOS confirmation and
prospective helper source pins remain required; scientific results are unchanged.

The [Python 3.13 fixture follow-up](results/PY313_FIXTURE_COMPATIBILITY_20261002.md)
fixes valid oversized-tar construction and iterable frame-key comparisons in
tests only, preserving consumer guards and exact allocation assertions.
All 241 focused local cases pass; new-patch Python 3.13 execution is pending.
The preceding actual macOS Python 3.13 log confirms all 135 cleanup/FAS/provider
cases, but retains 37 failures and does not establish complete CI or release.

The [Linux native diagnostic coverage](results/LINUX_NATIVE_DIAGNOSTIC_CI_20261002.md)
adds a capable CI runner for Linux-only constructor/worker checks, with a
fail-on-skip receipt gate. All 317 ordered local cases pass, including 20
native worker cases. New-runner confirmation is pending; synthetic parser
runtime binding does not establish native libc/parser replay or full release.
The preceding actual macOS Python 3.13 fast log confirms all 241 prior
archive/allocation/cleanup/FAS/provider cases; 28 other failures remain.
The first Linux job has 75 passes and one ignored-launcher failure. The
[explicit test staging follow-up](results/LINUX_NATIVE_LAUNCHER_STAGING_20261002.md)
retains that failure and gives 76 local passes; new-runner confirmation is pending.

The [frontier CPU replay correction](results/FRONTIER_SUM_REPLAY_20261002.md)
preserves historical aggregate addition order without relaxing exact replay.
All three archived method screens match under Python 3.12, including a
checkout-blocked copied child; 172 local tests pass. New source-bound recipes
must be frozen prospectively. This does not admit controlled timing or full release.

The [retained overhead metadata follow-up](results/RETAINED_OVERHEAD_METADATA_20261002.md)
binds only explicit test record paths, retaining hashes, sizes and complete
scientific payload equality. All 128 local cases pass; all 69 affected-module
cases also pass in a checkout-blocked copied tree. The preceding actual Python
3.13 log confirms the frontier correction, but has 14 other failures and does
not contain this new test-only fix. Full CI and release remain unfinished.

The [offline remote-path correction](results/REMOTE_NATIVE_PATHS_20261002.md)
preserves target paths without preparation-host symlink resolution, keeping
default local harness behavior and strict command equality. All 120 local
cases also pass in an original-checkout-blocked copied tree; the entire
15-field/27-run plan matches retained payloads. The preceding actual Python
3.13 log confirms all 69 metadata cases, but has ten other failures. New-patch
remote confirmation, source-inventory handoff and complete release remain open.

The [historical lineage replay follow-up](results/LINEAGE_FROZEN_REPLAY_20261002.md)
stages exact frozen dependencies for fresh-child tests instead of accepting
current helpers under old pins. All 134 focused cases pass; actual Python
3.12 reproduces three complete archived controls. Shared original Git objects
remain an explicit dependency, not a standalone archive. The preceding actual
Python 3.13 log confirms all 120 command/metadata cases, with eight other
failures. New-patch remote confirmation and complete release remain open.

The [corrected inventory test binding](results/CORRECTED_INVENTORY_BINDINGS_20261002.md)
retains source hashes and sizes while relocating four explicit test records.
All 122 focused local cases pass; 44 affected cases pass in a checkout-blocked
copy. Actual preceding Python 3.13 confirms all 21 frozen-lineage replay cases,
with six other failures. Raw exporter dependencies and full CI remain open;
no scientific scores change and no quiet window is required now.

The [explicit SwissTrees raw-source relocation workflow](results/SWISS_RAW_SOURCE_RELOCATION_20261002.md)
preserves every historical identity and raw check at new paths. All 187 focused
tests pass; both exporters produce byte-identical TSV/Markdown in a guarded
copied tree with the complete raw panels. This opt-in workflow does not rerun
native annotation extraction, provision CI data or authorize redistribution.

The [Linux host-reader CI placement](results/HOST_READER_PLATFORM_COVERAGE_20261002.md)
requires actual live execution on Linux, retaining portable parser/synthetic
coverage elsewhere. All 110 local cases pass with the live reader deliberately
not run; 34 copied cases pass with original and live-host reads blocked.
Actual Linux confirmation now has 111 passes, including the live reader;
see the [raw-regression handoff follow-up](results/SWISS_RAW_PYTEST_HANDOFF_20261002.md).
That follow-up provides explicit digest-pinned pytest options: 200 local cases
pass and all 27 affected cases pass in a copied tree, retaining default and
wrong-binding failures. Public CI data provisioning/rights remain unresolved;
neither result establishes controlled timing.

The [private raw archive restoration](results/SWISS_RAW_ARCHIVE_RESTORATION_20261002.md)
adds bounded, independently checksum-pinned input transfer. Both complete
archives restore under standard-library system Python 3.12; all 27 affected
raw regressions pass using those restored inputs with checkout reads blocked.
All 233 focused local cases pass. Raw archives stay private/local, with no
native annotation rerun, public upload, rights grant or full-release claim.

The [SwissTrees descriptive-table checker](results/SWISS_DESCRIPTIVE_TABLE_ARITHMETIC_20261002.md)
verifies 208 rows across all eight methods using exact rational arithmetic,
including 984 score/difference cells and explicit missingness. It runs with
standard-library Python and committed derived data; 49 local tests pass,
including an original-checkout-blocked copied child. It does not regenerate
raw annotations/scoring, bootstrap intervals or the existing review archive.

The [SwissTrees arithmetic archive](results/SWISS_DESCRIPTIVE_COMPONENT_20261002.md)
packages that unchanged checker and committed inputs. Its actual fresh
restoration under Python 3.12 reproduces all 208 rows and 984 cells with
original-checkout reads blocked; 76 combined local cases pass. The subsequent
actual macOS Python 3.10 fast log confirms all 76 checker/bundler cases pass,
but has 29 other failures. This derived-table component is not raw admission
or the full study release.

The [provider-checksum follow-up](results/PROVIDER_CHECKSUM_PORTABILITY_20261002.md)
preserves BSD checksum semantics while removing a macOS-incompatible flag.
Malformed output and execution failures cannot become verified matches;
138 focused local tests pass. Historical source verification is not rerun or
repinned, and data rights, controlled timing and full release remain open.
Subsequent actual macOS Python 3.10 fast CI confirms all 138 cases, but still
has 29 other failures, including an owned-timeout EPERM failure.

The [native FAS source reference](results/FAS_NATIVE_SOURCE_REFERENCE_20261002.md)
retains the exact upstream scorer and license for explicit AST comparison tests.
All 33 FAS cases pass in a copied tree with original-checkout reads blocked;
65 combined local cases pass. No source fallback, biological data substitution,
full FAS execution or change to historical source admission follows.
The subsequent actual macOS Python 3.13 fast log confirms all 65 FAS/provider
cases pass, but has 37 other failures. This is not full remote-suite success.

The FastOMA and OrthoMCL launch examples below describe historical workflows.
They are **not** the corrected-input publication execution recipe. Do not
launch them to resume the active corrected-release jobs or transfer their
historical scores to the corrected input release. Consult the
[progress ledger](results/PUBLICATION_PROGRESS.md) and verify live scheduler
state before any submission. An observation timeout is not a failed job.

## FastOMA comparator

Run FastOMA 0.3.5 on the complete 78-proteome QfO 2020 input and chain the
result into the retained six-metric assessment with:

```bash
benchmark_tools/submit_fastoma_qfo.sh
```

The inference runner uses the cached LUCA OMAmer database, emits bare UniProt
accessions, and forces FastOMA's pairwise ortholog output for the 78-species
dataset. FastOMA requires a supplied species tree, so this comparator uses the
OrthoFinder 3.1.5 tree inferred from the same input proteomes. This tree does
not contain QfO benchmark labels. The runner records the tree checksum and
refuses to overwrite an existing production result. OMAmer searches start
with 24 GB each and are capped at 28 concurrent tasks to stay within the
default 700 GB job allocation.

If FastOMA or its executor fails after producing reusable Nextflow work,
resume it and recreate the dependency-chained scoring job with:

```bash
RESUME=true benchmark_tools/submit_fastoma_qfo.sh
```

FastOMA 0.3.5 rejects a lowercase `uniprot` identifier transformer during
final collection even though all inference tasks have finished. For a legacy
run that reached this exact failure, recover its staged 271-batch collection
without rerunning inference:

```bash
COLLECT_WORK_DIR=/absolute/path/to/failed/collect_subhogs/workdir \
  benchmark_tools/submit_fastoma_qfo_recovery.sh
```

The completed FastOMA result and the current cross-tool QfO/OrthoBench tables
are recorded in `results/EXTERNAL_TOOL_COMPARISON_20260904.md` and its paired
machine-readable JSON file.

## OrthoMCL comparator

Run the restartable OrthoMCL 1.4 QfO benchmark and queue official scoring:

```bash
benchmark_tools/submit_orthomcl_qfo.sh
```

This retains OrthoMCL 1.4's legacy BLAST 2.2.13 search and MCL 02-063
clustering defaults. The runner replaces the impractical serial BPO parser
with the byte-compatible streaming converter validated against native output,
and uses checkpointed process parallelism for independent species pairs. A
rerun with the same output directory resumes at the first incomplete stage.

## Production OrthoHMM performance harness

Use the production CLI benchmark harness for runtime and memory measurements:

```bash
python benchmark_tools/benchmark_production.py \
  bacterial_scaling/subsets/n005 \
  /tmp/orthohmm_n005 \
  /tmp/orthohmm_n005.metrics.json \
  --cpu 32 \
  --accuracy-profile high_sensitivity
```

The JSON records exact commands, input/output checksums, stage wall time,
sampled process-tree RSS, search candidates, significant hits, graph edges,
and orthogroup counts. The bacterial SLURM runner uses this harness rather
than maintaining a second OrthoHMM implementation.

The harness defaults to `--accuracy-profile standard`. Select
`high_sensitivity` explicitly when measuring the bounded k=4 search,
multi-pass graph inference, and strict profile expansion workflow.

## OrthoBench stage diagnostics

Use the offline diagnostic to measure whether reference relationships survive
search hits, graph construction, clustering, and refinement:

```bash
python benchmark_tools/orthobench_stage_diagnostics.py \
  --refogs /path/to/Open_Orthobench/BENCHMARKS/RefOGs \
  --hits-pickle mc100_hits=benchmarks/results/hits_BLOSUM62_mc100.pkl \
  --pairs production_edges=/path/to/orthohmm_edges.txt \
  --clusters production_final=/path/to/orthohmm_orthogroups.txt \
  --official-benchmark /path/to/Open_Orthobench/BENCHMARKS/benchmark.py \
  --json /path/to/stage_diagnostics.json
```

The tool records input hashes and reports both direct reference-pair recovery
and recovery through connected components. Reference labels are read only by
this post-run process; they are never passed into OrthoHMM inference. Pickle
input is provided only for trusted local historical caches.

## High-sensitivity graph replay

Validate production graph and refinement code against a trusted normalized-hit
checkpoint without rerunning the all-to-all sequence search:

```bash
python benchmark_tools/replay_high_sensitivity.py \
  --hits-pickle benchmarks/results/hits_BLOSUM62_mc100.pkl \
  --fasta-directory /path/to/Open_Orthobench/BENCHMARKS/Input \
  --output-directory /tmp/orthohmm_high_replay \
  --official-benchmark /path/to/Open_Orthobench/BENCHMARKS/benchmark.py \
  --json /tmp/orthohmm_high_replay.json
```

The replay uses no reference labels during inference. Its JSON records source
and input checksums, graph counts, runtime, peak process RSS, output checksums,
and optional post-inference OrthoBench scores. Omit `--fasta-directory` to
validate only RBNH, singleton reassignment, and production refinement.

Use `--profile-iterations 2` only to evaluate iterative profile rebuilding.
The production-equivalent default remains one pass. The recorded two-pass
OrthoBench experiment doubled profile-search work and reduced refined F-score
from 70.4 to 70.2, so iterative refinement is not a production recommendation.

Use `python -m benchmark_tools.selective_jackknife_holdout` to test profile
threshold calibration in the presence of mixed, duplicated-species profiles.
The replay flag `--jackknife-single-copy-profiles` applies the jackknife only
when every profile member comes from a distinct species. This mode is separate
from the rejected global `--jackknife-profile-thresholds` ablation.

## Phylogeny-stage replay

Replay the optional phylogeny stage from a frozen candidate-cluster file:

```bash
python benchmark_tools/replay_phylogeny.py \
  --fasta-directory /path/to/proteomes \
  --candidate-clusters /path/to/orthohmm_edges_clustered.txt \
  --output-directory /tmp/orthohmm_phylogeny_replay \
  --json /tmp/orthohmm_phylogeny_replay.json \
  --species-tree-mode infer --cpu 32 \
  --root-rule supported_children --pair-rule positive_paralogy \
  --official-benchmark /path/to/OrthoBench/benchmark.py
```

The JSON records source/input checksums, process-tree memory, stage runtime,
reconciliation counts, output checksums, and optional post-run scores. This is
a development harness; final claims still require fresh end-to-end runs.
For a reconciliation-only comparison, pass `--checkpoint-source` with a prior
replay output. The harness hard-links immutable raw-tree checkpoints when the
filesystem permits, recomputes rooting and reconciliation in a distinct output
directory, and records the source path and checkpoint-hit counts.

## Permissive HMM candidate graphs

Audit a frozen candidate partition before spending time on gene-tree
inference. The diagnostic verifies that candidates only combine immutable seed
families and reports the reference-pair recall ceiling, reference contamination,
family growth, and per-RefOG component recovery:

```bash
python benchmark_tools/candidate_family_diagnostics.py \
  --seed-clusters /path/to/seed_clusters.txt \
  --candidate-clusters /path/to/candidate_clusters.txt \
  --refogs /path/to/OrthoBench/RefOGs \
  --json /tmp/candidate_audit.json
```

Build candidate superfamilies from a trusted normalized-hit checkpoint before
testing phylogenetic reconciliation:

```bash
python benchmark_tools/build_permissive_candidates.py \
  --hits-pickle benchmarks/results/hits_BLOSUM62_mc100.pkl \
  --output-directory /tmp/orthohmm_candidates \
  --json /tmp/orthohmm_candidates/result.json \
  --threshold-factor 0.8 --inflation 1.2 \
  --mcl /path/to/mcl --cpu 32 \
  --official-benchmark /path/to/Open_Orthobench/BENCHMARKS/benchmark.py
```

The support factor scales each gene's weakest reciprocal-best normalized-hit
threshold. Graph construction is label-blind, MCL receives binary edge weights,
and omitted isolates are restored before the candidate file is finalized. The
JSON records parameters, checksums, stage runtimes, peak process-tree memory,
graph/family counts, and any post-inference OrthoBench score.

To preserve OrthoHMM's precise seed groups and merge only mutually
best-supported group pairs, run:

```bash
python benchmark_tools/merge_candidate_superfamilies.py \
  --hits-pickle benchmarks/results/hits_BLOSUM62_mc100.pkl \
  --seed-clusters /path/to/orthohmm_edges_clustered.txt \
  --output-directory /tmp/orthohmm_superfamilies \
  --json /tmp/orthohmm_superfamilies/result.json \
  --max-family-genes 500 --iterations 2 \
  --official-benchmark /path/to/Open_Orthobench/BENCHMARKS/benchmark.py
```

This alternative aggregates directed HMM evidence between seed groups and
joins only reciprocal-best partners. It is intended to avoid giant components
while recovering split remote homologs for subsequent reconciliation. Bounded
iterations let a newly merged group recruit another mutually best-supported
fragment without relaxing the evidence rule.

The experimental `--method satellite` mode relaxes mutual-best ranking only
for a bounded small source group. It still requires bidirectional HMM support,
a source-side best-target margin, size and component caps, and an optional
species-overlap limit. This is intended for fragments whose larger anchor has
another stronger bridge:

```bash
python benchmark_tools/merge_candidate_superfamilies.py \
  --method satellite --hits-pickle /path/to/hits.pkl \
  --seed-clusters /path/to/orthohmm_edges_clustered.txt \
  --output-directory /tmp/orthohmm_satellite_candidates \
  --json /tmp/orthohmm_satellite_candidates.json \
  --max-satellite-genes 5 --max-satellite-ratio 0.75 \
  --min-margin 1.5 --max-satellites-per-anchor 2 \
  --min-coverage 0.5 --min-normalized-support 0.02
```

Validate that bounded behavior independently with:

```bash
python benchmark_tools/candidate_superfamily_holdout.py \
  benchmark_tools/results/candidate_superfamily_holdout.json
```

The synthetic graph contains two-to-five-fragment true families and a cycle of
weaker asymmetric cross-family decoys. It reports pair precision and recall for
one pass and four bounded passes across fixed holdout seeds.

Root-HOG duplication evidence can be replayed without rebuilding alignments or
trees:

```bash
python benchmark_tools/replay_root_hog_rules.py \
  --fasta-directory /path/to/proteomes \
  --candidate-clusters /path/to/candidates.txt \
  --phylogeny-output /path/to/run/orthohmm_phylogeny \
  --output-directory /tmp/root_hog_rules \
  --json /tmp/root_hog_rules/result.json \
  --official-benchmark /path/to/Open_Orthobench/BENCHMARKS/benchmark.py
```

The production default remains `supported_children`; the replay compares that
rule with progressively broader root-level duplication evidence from the same
rooted trees.

Before tuning candidate-family or reconciliation parameters, freeze the
OrthoBench development and validation RefOGs without reading their members:

```bash
python benchmark_tools/create_orthobench_partition.py \
  --refogs /path/to/Open_Orthobench/BENCHMARKS/RefOGs \
  --json benchmark_tools/results/orthobench_partition.json
```

Score a prediction on only one frozen partition with the official
equal-RefOG pairwise formula:

```bash
python benchmark_tools/score_orthobench_partition.py \
  --predictions /path/to/root_hogs.txt \
  --refogs /path/to/Open_Orthobench/BENCHMARKS/RefOGs \
  --partition-json benchmark_tools/results/orthobench_partition.json \
  --partition development --json /tmp/development_score.json
```

Synthetic phylogeny holdouts accept independent reconciliation policies. For
example, evaluate the policy that preserves pairs at mapping-only conflicts
while retaining the conservative root-HOG rule with:

```bash
python benchmark_tools/phylogeny_holdout.py /tmp/phylogeny_holdout.json \
  --work-directory /tmp/phylogeny_holdout \
  --aligner /path/to/mafft --tree-builder /path/to/FastTree --cpu 4 \
  --root-rule supported_children --pair-rule positive_paralogy
```

After predictions are frozen, stratify reconciliation changes by RefOG copy
number and alignment-derived sequence identity:

```bash
python benchmark_tools/analyze_phylogeny_changes.py \
  --refogs /path/to/Open_Orthobench/BENCHMARKS/RefOGs \
  --fasta-directory /path/to/Open_Orthobench/BENCHMARKS/Input \
  --baseline-clusters /path/to/orthohmm_edges_clustered.txt \
  --phylogeny-root-hogs /path/to/orthohmm_root_hogs.tsv \
  --alignment-directory /tmp/refog_alignments \
  --aligner /path/to/mafft --cpu 4 \
  --json /tmp/phylogeny_changes.json
```

This diagnostic reads benchmark labels only after inference. It verifies that
root HOGs remain subsets of their source families, distinguishes splitting
from cross-family merging, counts retained and removed reference relations,
and reports copy-number and median pairwise-identity strata.

## OrthoFinder comparator

This directory pins the OrthoFinder comparator to version 3.1.5 and keeps new
results separate from the historical 2.5.5 outputs.

## Installation

```bash
benchmark_tools/install_orthofinder.sh
```

The installer verifies the upstream archive checksum and the installed version.
Set `SOFTWARE_ROOT` or pass an explicit installation directory when the shared
workspace layout is unavailable.

## Benchmarks

Submit the serialized OrthoBench and QfO workflow with:

```bash
benchmark_tools/submit_orthofinder_v3_benchmarks.sh
```

The QfO runner writes both raw native OrthoFinder pairs (`pairs.tsv`) and pairs
whose IDs occur in the QfO mapping (`pairs.qfo.tsv`). Only the filtered file is
submitted to the official assessment pipeline. The scoring workflow retains
VGNC, SwissTrees, TreeFam-A, EC, GO, and FAS.

To isolate the contribution of OrthoFinder's phylogenetic stages, convert the
pre-phylogenetic DIAMOND/MCL checkpoint retained in `WorkingDirectory`:

```bash
python benchmark_tools/orthofinder_mcl_to_orthogroups.py \
  WorkingDirectory/clusters_OrthoFinder_I1.2.txt_id_pairs.txt \
  WorkingDirectory/SequenceIDs.txt \
  orthogroups_sequence_only.txt
```

This output precedes MSA, gene-tree, species-tree, reconciliation, and HOG
delineation. For QfO it is converted to cross-species co-cluster pairs with
`qfo_benchmark/og_to_pairwise.py`, the same representation used for every
group-only method.

Useful environment overrides include `ORTHOFINDER_BIN`, `SOFTWARE_ROOT`,
`ORTHOBENCH_ROOT`, `QFO_PIPELINE`, `QFO_MAPPING`, `NEXTFLOW_BIN`, and
`JAVA_HOME_17`. Set `RESUME=1` when retrying an interrupted QfO assessment so
Nextflow reuses completed tasks instead of deleting its work directory.

Summarize completed QfO assessments with:

```bash
python benchmark_tools/qfo_summarize_scores.py \
  qfo_benchmark/scoring/orthofinder_v3_diamond
```

Run the production OrthoHMM high-sensitivity QfO workflow separately with:

```bash
benchmark_tools/submit_orthohmm_high_sensitivity_qfo.sh
```

The inference job uses `benchmark_production.py`, records production stage
metrics and manifests, converts final orthogroups to cross-species QfO pairs,
and filters identifiers against the official mapping. Official six-metric QfO
scoring starts only after successful inference.

After freezing a phylogeny-aware candidate, submit its QfO workflow with:

```bash
benchmark_tools/submit_orthohmm_phylogeny_qfo.sh
```

This run freezes `high_sensitivity` plus the `satellite_v2` HMM candidate
profile, minimum-variance internal species-tree rooting, the `species_overlap`
root-HOG rule, and the `positive_paralogy` pair rule. It infers its species tree
internally and submits OrthoHMM's native reconciled cross-species pairs to QfO.
It does not expand root HOGs into cliques.

If a fresh production run reaches a verified candidate-cluster checkpoint but
fails during phylogeny, `submit_resume_orthohmm_phylogeny_qfo.sh` resumes only
that stage. The restart harness verifies and preserves the failed-run
checkpoint, retains both metrics files, and emits `result.json` with combined
wall time and peak process-tree memory before submitting QfO scoring.

Run the secondary Three Kingdoms validation with:

```bash
benchmark_tools/submit_orthohmm_high_sensitivity_three_kingdoms.sh
```

This uses the same production benchmark harness and accuracy profile, then
scores the final orthogroups against the local BUSCO-derived reference. The
workflow expects exactly 12 FASTA files and refuses to replace an existing
result directory.

## Bounded Copy-Branch Hill-Climb

`hillclimb_copy_branch_ablation_20261013.json` prespecifies separate deferral of
the large-copy and medium-panel copy-splitting branches, retaining the original
500-gene ceiling and all eligibility safeguards. The experimental profiles are
`phylogeny_deferred_large` and `phylogeny_deferred_medium`; neither changes
production defaults. HMM-derived graph inputs, satellite expansion, species-
overlap evidence, and positive-paralogy inference remain matched to the control.

Prepare each dataset in a new namespace using an immutable source checkout:

```bash
CUDA_VISIBLE_DEVICES='' OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=16 PYTHONHASHSEED=0 \
python "$SOURCE/benchmark_tools/run_copy_branch_ablation.py" prepare \
  --repo "$REPO" --source-root "$SOURCE" --output "$OUTPUT" \
  --dataset "$DATASET" --parent-root "$PARENT" --parent-job "$PARENT_JOB"
```

For QfO, additionally pass `--fixed-root` pointing to the admitted fixed-tree
copy-deferral graph experiment. Submit the frozen Slurm launcher with source
root and execution-plan path as its two arguments, after the parent terminates.
The worker independently admits the parent using scheduler accounting and
content readback before constructing candidates. The recorded launch receipt
is `results/copy_branch_ablation_launch_20261013.json`.

The large-copy branch cannot activate on the 12-species OrthoBench and Three
Kingdoms panels. Their large-only arm must exactly reproduce control candidates
and constraints; the medium-only arm must reproduce combined deferral. After
downstream-source equivalence checks, these arms reuse the corresponding native
inference and independently re-read scores. They are structural controls, not
new accuracy gains or inference timings. QfO's two new native arms run
sequentially with candidate-specific inferred species trees.

Add `--copy-branch-root` to `report_hillclimb.py` to include the experiment's
`OrthoBench`, `QfO`, and `ThreeKingdoms` subdirectories. The independent reader
verifies parent admission, exact reuse, native gene/pair coverage, scoring
handoffs, and executed deferral traces. All eight endpoints remain separate;
the 0.005 decline threshold is a development screen, not an equivalence claim.
Missing or failed cells stay scoreless. Local shared-host timings and cached
HMM graph replays must not be ranked as isolated fresh full-pipeline runtimes.

For a completed own-inferred deferral round, localize its retained GO/EC scored
pair changes against native reconciliation-confidence tables with:

```bash
python benchmark_tools/diagnose_hillclimb_pair_confidence.py \
  --receipt "$INDEPENDENT_READBACK/summary.json" --output "$NEW_DIAGNOSTIC"
```

This checks receipt-bound native coverage and raw endpoint counts/means without
rerunning inference or scoring. High/medium event labels are not calibrated
correctness probabilities; annotated subsets and shared-pair conditioning do
not replace the original benchmark endpoints.

## Bounded Large-Copy Components

`hillclimb_large_copy_components_20261013.json` prespecifies the opt-in
`phylogeny_components` refinement profile. It retains within-family components
at the existing edge-weight threshold of 1.5 only when each component has at
most 500 genes. Oversized components and isolated genes retain singleton
fallback. The unchanged large-copy eligibility predicate controls activation;
medium-panel refinement is unchanged. Components supply candidate context,
not ortholog assignments.

The native CLI/API requires fresh built-in high-sensitivity search, Leiden,
phylogenetic reconciliation, `satellite_v2` candidates, `species_overlap`
rooting, `positive_paralogy` pair inference, `high_confidence_pair` satellite
membership and an inferred `min_variance` species tree. Whole-family deferral
cannot be combined with component retention. Production defaults are unchanged.

The all-family feasibility probe and independent union-find readback are
retained in `results/large_copy_component_feasibility_v3_20261013.json` and
`results/large_copy_component_independent_readback_20261013.json`. They do not
measure accuracy. Finish and independently admit the existing copy-branch
vectors before launching the next component inference round; then freeze
source, inputs, scoring and resource settings and evaluate all eight endpoints.

After all four branch vectors are independently admitted, prepare one fresh
component namespace per dataset with the frozen runner:

```bash
CUDA_VISIBLE_DEVICES='' OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=16 PYTHONHASHSEED=0 \
python "$SOURCE/benchmark_tools/run_large_copy_components.py" prepare \
  --repo "$REPO" --source-root "$SOURCE" --output "$OUTPUT" \
  --dataset "$DATASET" --branch-root "$BRANCH_ROOT" \
  --comparison "$COMPLETE_BRANCH_COMPARISON"
```

Preparation requires the actual successful terminal branch job and verified
cells for all eight endpoints; a live or partial report cannot trigger new
inference. Submit `run_large_copy_components.slurm` with source root and
execution-plan path. The worker reproduces the original default seed,
candidates and constraints before evaluating the new profile. QfO uses
candidate-specific inferred trees; inactive narrow panels reuse native results
only after exact candidate checks and reread their original scores. Identical
QfO candidate/constraint bytes also avoid duplicate native inference.

Independently reconstruct candidates and reread terminal outputs with:

```bash
CUDA_VISIBLE_DEVICES='' OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=16 PYTHONHASHSEED=0 \
python "$SOURCE/benchmark_tools/readback_large_copy_components.py" \
  --execution-plan "$OUTPUT/execution_plan.json" --output "$NEW_READBACK"
```

Candidate reconstruction is a separate cached validation cost, not another
native inference run. Keep its cost separate from shared-host native/scoring
timings, and do not interpret an inactive panel's unchanged score as evidence
that the large-copy branch generalizes.

A separately frozen repaired reader may pass `--reader-source-root "$READER"`
when reading an older execution plan. The reader verifies that every scientific
Python module and benchmark helper except the reader itself is byte-identical
to the producer checkout; the producer's input and source checks remain intact.
