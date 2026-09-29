# Publication Reproduction Guide

The [eight-method OrthoBench uncertainty extension](results/OB_COMPLETE_UNCERTAINTY_RESULT_20260928.md)
reports all 21 paired contrasts under a frozen exploratory protocol. It
preserves historical intervals separately and does not establish independent
confirmation or validate the family-exchangeability assumption.

The [candidate wheel dependency audit](results/RELEASE_WHEEL_DEPENDENCIES_20260928.md)
checks declared runtime dependencies against both hash-pinned wheel sets on
the recorded interpreter. It does not establish native/OS closure or security.

Status: 27 September 2026, incomplete working package. This guide routes
readers to executed workflows and their evidence; it is not an archival
release or an assertion that every requirement is complete. Direct commands
below recompute statistics or audit retained local artifacts; they do not
submit scheduler jobs. Linked native batch recipes have separate execution
requirements and must not be confused with these audit commands.

## Execution Status

These are distinct executed paths, not interchangeable evidence:

| Path | Verified scope | Does not establish |
| --- | --- | --- |
| [Full native recovery, jobs 22326/22327](results/FULL_RECOVERY_ORTHOBENCH_RESULT_22326.md) | Fresh search/trees; all OrthoBench groups and 70 family scores reproduced | Independent accuracy or controlled timing |
| [Restored scoring archive](results/OB_SCORING_ARCHIVE_20260927.md) | Full OrthoBench scoring outside the checkout, using supplied raw inputs and reader | Fresh acquisition, installation or inference |
| [Integrated fixture](results/INTEGRATED_WORKFLOW_20260927.md) | Eight stages on 16 genes, with fresh inference/reader environments | Full-data validation |
| [Reconstructed-base fixture](results/RECONSTRUCTED_BASE_FIXTURE_20260927.md) | Same fixture using a separately acquired base Python; no original Python/repository prefix in its execution trace | Cross-host restoration, OS isolation or complete security/rights clearance |
| [Full integrated job 22337](results/INTEGRATED_FULL_OB_RESULT_22337.md) | Eight stages completed; separate admission reproduced all groups and 70 family scores | Controlled timing, independent accuracy or cross-host restoration |
| [Reconstructed-base full job 22376](results/RECONSTRUCTED_FULL_OB_PROTOCOL_20260929.md) | Submitted under a frozen protocol; live installation check passes for 11 inference and five reader packages | Native completion, scientific admission or reproduction equality are not yet established |

The reconstructed-base full attempt has separate execution and scientific
checkers, `benchmark_tools.verify_reconstructed_full_ob_execution` and
`benchmark_tools.admit_reconstructed_full_ob`. They are bound to job 22376
and its exact plan; do not substitute another job ID. See the newest entries
in the [progress ledger](results/PUBLICATION_PROGRESS.md) and query the
scheduler for current state. Its installation-only payload audit does not
replace terminal admission. The comparison baseline is admitted job 22337;
root groups are compared without labels and native pair/event/hierarchy
files are compared byte-for-byte, retaining any differences. No outcome of
this shared-host correctness run qualifies as controlled timing.

For a new execution, acquire and verify the selected inputs and local assets
first, then follow the integrated controller's documented command with fresh
output paths. Reference files remain separate from inference FASTA inputs.
The reconstructed-base receipt records a specific 19-package Conda selection
plus a hash-pinned pip wheel overlay; it is not a solver-complete Conda lock
or a recommendation to clone the development environment. Do not substitute
an unvalidated dependency update into a frozen scientific run.

Job 22337's admission command is deliberately bound to its exact job and
plan. It must not be used to approve a different execution by changing the
job identifier or treating the fixture's success as dataset-scale evidence.
Consult current scheduler state before calling it; the protocol's historical
RUNNING observation does not establish current liveness or completion.

## Method And Results

- [Corrected GO/EC scored-pair panel](results/QFO_SCORED_PAIR_PANEL_20260927.md)
  provides the executable eight-method batch command and all 56 pairwise
  comparisons. Shared pairs have identical serialized scores; original means
  retain each method's own eligible-pair set and denominator. This is not a
  family-aware confidence interval or a replacement benchmark endpoint.

- [Completed QfO ordering comparison](results/QFO_CANONICAL_RESULT_22333.md)
  distinguishes exact fresh retained-order reproduction from changed canonical
  predictions: 110 pairs lost, 85 gained and seven shared-pair annotation
  changes, with identical species trees. The [independently admitted scoring](results/QFO_CANONICAL_ASSESSMENT_RESULT_22336.md)
  finds five unchanged endpoint scores and a small, sampling-confounded FAS
  difference. Both evaluated rows are preserved; no default was selected by
  scores and no accuracy result was transferred from another dataset.

- [Fresh reader-only environment](results/FRESH_READER_RUNTIME_20260927.md)
  runs the [exported independent validators](results/RELOCATED_INDEPENDENT_READERS_20260927.md)
  outside the checkout using five locally installed, hash-locked packages.
  All four fixture readers pass, with 2,700 package files checked before/after.
  The subsequent [Biopython security upgrade](results/READER_SECURITY_UPGRADE_20260927.md)
  validates the separate 1.87 lock with unchanged scientific fields and 2,702
  audited package files. Use that patched reader lock for new installations;
  retain the earlier lock only as historical provenance.
  This is separate from the inference runtime; base Python/OS remain shared.

- [Full OrthoBench scoring outside the checkout](results/PORTABLE_OB_SCORE_20260927.md)
  connects separately acquired raw inputs, exported readers and the patched
  reader environment. All 70 reference-family score records reproduce exactly
  from the retained full prediction partition. This is full-data scoring,
  not a new native run or a complete portable acquisition-to-inference workflow.
  Its [deterministic scoring archive](results/OB_SCORING_ARCHIVE_20260927.md)
  has also been restored and executed outside the checkout, with byte-identical
  full score output. Raw upstream inputs and the reader interpreter must be
  supplied separately; the local archive is not a public release.

- [Integrated installation-to-scoring controller](results/INTEGRATED_WORKFLOW_20260927.md)
  executes all eight stages on the installation fixture with separate fresh
  environments. It also detects and corrects two historical absolute MAFFT
  convenience links in a new asset copy. Full-OrthoBench inputs are pinned,
  and [full-data job 22337](results/INTEGRATED_FULL_OB_RESULT_22337.md)
  completed under a protocol frozen before submission. Independent admission
  reproduced all 59,770 groups and 70 family scores, with byte-identical root
  partitions. This run used the original base Python on the shared host.
  The [independent post-completion admission command](results/INTEGRATED_FULL_OB_ADMISSION_20260927.md)
  is implemented and tested; it refuses a still-running or failed job.
  The [exact integrated-wheel notice texts](results/INTEGRATED_DEPENDENCY_NOTICES_20260927.md)
  have been collected and verified separately for release review, not redistribution clearance.

- [Latest full unit regression](results/PUBLICATION_TEST_REFRESH_20260929.md)
  records 11,885 passes and 12 skips at revision `cc3df68a`. Ten skips are
  opt-in native checks; two are Selectome module-collection skips for absent
  `sqlglot`. A separate run with the retained parser dependency passes all
  13 tests in those two modules; this is not a combined full-suite result.
  The [10,874-pass run](results/PUBLICATION_TEST_REFRESH_20260927_v3.md)
  remains historical evidence at `3eb446e`.
  The 10,819-pass run remains historical evidence at `75d42e7`.
  The earlier 10,729-pass run remains historical evidence at `9914060`.
  The [separate installed-native follow-up](results/PUBLICATION_NATIVE_TEST_REFRESH_20260927.md)
  passed the earlier run's ten skipped tests; it is not automatically a
  native-suite pass at every subsequent revision. The later
  [ELF inventory](results/INTEGRATED_WHEEL_ELF_20260927.md) has 41 focused
  passing focused tests and is included in the full suites, alongside
  base-archive staging and GO/EC scored-pair tools. Preserve these scopes
  and the 22 frozen-source warnings; tests do not establish biological
  accuracy or controlled performance.

- [Separate base-runtime reconstruction](results/RECONSTRUCTED_BASE_FIXTURE_20260927.md)
  supplies fresh Python for installation and execution of the integrated
  fixture. Its 6,229 audited bootstrap/inference/reader wheel payloads match.
  The [Conda payload audit](results/RECONSTRUCTED_CONDA_PAYLOADS_20260927.md)
  covers 6,664 declared entries; [archive-based replay](results/RECONSTRUCTED_PREFIX_REPLAY_20260927.md)
  reproduces all 195 prefix-rewritten files. The one Python bytecode-cache
  serialization difference remains recorded with matching compiled code
  fields, not silently recategorized as archive-byte equality. These checks
  do not replace OS provenance or validate unobserved runtime paths.
  The [archive-staging command](results/BASE_ARCHIVE_STAGING_20260927.md)
  recreates the checked local explicit file from a supplied archive cache;
  22 focused tests and two real staging executions pass. It performs no
  download or installation and postdates the latest full unit suite.

- [Observed search libraries](results/INTEGRATED_SEARCH_LIBRARIES_20260927.md)
  and [base package attribution](results/INTEGRATED_BASE_RUNTIME_PACKAGES_20260927.md)
  distinguish wheel bytes from base Python/system providers. The static ELF
  inventory lists declared dependencies, not runtime resolution; the live
  snapshot covers two search-stage processes, not every later stage. Both
  [wheel notices](results/INTEGRATED_DEPENDENCY_NOTICES_20260927.md) and
  [base-package notices](results/RECONSTRUCTED_BASE_NOTICES_20260927.md) have
  verified local exports. Component/source obligations and redistribution
  clearance remain open. The base mutex package's missing declaration is
  now scoped to a [metadata-only archive](results/BASE_MUTEX_PAYLOAD_20260929.md):
  its exact payload and installed-file lists are empty. This finding is not
  redistribution clearance for the archive or the actual GCC runtime.

- [Canonical assessment source archive](results/QFO_CANONICAL_SOURCE_ARCHIVE_20260927.md)
  preserves all 762 distinct repository-contained Python files pinned by the
  completed scoring plan. Later source edits must not be silently accepted
  by historical validators. The archive is local and source-only, not a
  complete scoring runtime or deposited release.

- [Fresh relocated runtime](results/FRESH_RELOCATED_RUNTIME_20260927.md)
  validates a new offline-installed environment and relocated MAFFT/FastTree
  on the small full-pipeline fixture, with package-byte audits and native
  file-access tracing. This still shares the host's base Python/system
  libraries and is not full-data or cross-platform validation.

- [Relocated recovery entrypoint](results/RELOCATED_ENTRYPOINT_FIXTURE_20260927.md)
  exports the exact committed minimal runtime harness and verifies a 16-gene
  execution outside the checkout with four scientific readers. It reuses the
  existing Python environment and external tools, so it is not a portable
  installation or full-data reproduction on another host.

- [Completed full OrthoBench recovery](results/FULL_RECOVERY_ORTHOBENCH_RESULT_22326.md)
  reproduces all 59,770 root groups and all 70 family-score records from raw
  FASTA in the new recovery installation. Native 22326 and independent audit
  22327 completed successfully. This is same-host reproducibility, not new
  accuracy evidence or controlled comparative timing. Manuscript review v25
  includes this completion, later QfO ordering results and bounded
  installation/scoring evidence; subsequent base-runtime checks are linked here.

- [Recovery advisory range review](results/RECOVERY_ADVISORY_REVIEW_20260927.md)
  matches the live recovery package inventory to its install report and checks
  the 11 fresh repository-alert ranges. None affects those installed versions;
  this is not comprehensive security clearance or closure of repository alerts.

- [Three Kingdoms OrthoMCL retained-input trace](results/THREE_KINGDOMS_ORTHOMCL_INPUTS_20260927.md)
  verifies all twelve native input copies, all 443,217 merged search sequences
  and the genome map against staged inputs. It does not establish immutable
  historical consumption, complete BLAST/BPO provenance or genome-wide accuracy.

- [Three Kingdoms fixed-release reacquisition](results/THREE_KINGDOMS_PROVIDER_REACQUISITION_20260928.md)
  freshly downloads all eleven fixed Ensembl sources and verifies exact size
  and SHA256 equality to retained compressed inputs. The moving Xenopus URL
  remains unresolved; this is not redistribution clearance or per-tool provenance.

- [Full OrthoBench recovery protocol](results/FULL_RECOVERY_ORTHOBENCH_PROTOCOL_20260926.md)
  freezes the 251,378-gene end-to-end run through the new explicit entrypoint.
  Job 22326 was submitted after protocol commit; its completed evidence is
  linked above. Original submission and live-state receipts remain historical.

- [Explicit full-pipeline fixture](results/PUBLICATION_FULL_PIPELINE_FIXTURE_20260926.md)
  runs raw FASTA through the frozen installed HMM pipeline with candidate-only
  canonical ordering and fresh inferred phylogeny. Four independent readers
  validate the 16-gene result. The separately executed full-data validation is
  linked above; its score was recomputed rather than transferred from a fixture.

- [Fresh recovery-install phylogeny fixture](results/PUBLICATION_RECOVERY_PHYLOGENY_20260926.md)
  validates inferred species/gene trees and reconciliation in the new hash-locked
  venv with four independent scientific readbacks. This 16-gene same-host test
  does not transfer the full OrthoBench recovery score to the new installation
  or complete portable full-pipeline reproduction.

- [Matched-recall simulation graph control](results/MATCHED_GRAPH_RESULT_20260926.md)
  separates search calibration from downstream orthology scoring. All 70 search
  datasets and 70 graph arms completed; reporting covers 35 datasets in seven
  conditions. It supports a bounded initial-search result, not full-pipeline
  superiority or matched real-data sensitivity. To recompute from the retained
  admitted local artifacts, use fresh output paths:

  ```bash
  python -m benchmark_tools.score_matched_graph \
    --readback benchmark_tools/results/matched_graph_readback_v2_20260926.json \
    --output /absolute/fresh/matched_graph_scores
  python -m benchmark_tools.plot_matched_graph \
    --results /absolute/fresh/matched_graph_scores/results.json \
    --output /absolute/fresh/matched_graph_figure
  ```

  This rechecks pinned local truth, input and prediction files; it does not
  acquire missing datasets, rerun inference or establish portable native
  reproduction. Exact native commands and resources remain in the
  [search submission](results/search_sensitivity_submission_20260926.json) and
  [graph submission](results/matched_graph_submission_20260926.json).
  The v2 readback adds exact environment/artifact bindings with identical cell
  evidence. The historical v1 report is retained; its auditor source hash
  predates the current strengthened auditor, so use v2 in the current checkout.

- [Current cross-dataset score table](results/CURRENT_BENCHMARK_SCORES_20260926.md)
  combines retained OrthoBench, corrected QfO and supplementary Three Kingdoms
  scores without mixing QfO input releases or averaging across datasets.
  Version 2 adds pinned five-comparator OrthoBench readback/upstream evidence
  without changing the score TSV or claiming historical native provenance.

- [Refreshed extended manuscript review](results/MANUSCRIPT_RENDER_REVIEW_20260928_v27.md)
  provides dated HTML/PDF and a checked local-asset inventory. It remains
  a working preview, not a standalone archive or completed journal layout.
  Version 27 adds the complete descriptive OrthoBench figure and refreshes
  the intervening source text; selected new pages were visually inspected,
  not the complete document. Earlier snapshots remain retained.
  None adds native VGNC confidence intervals or controlled timing.
- [Corrected six-endpoint comparison figure](results/CORRECTED_QFO_ENDPOINT_FIGURE_20260926.md)
  retains all eight methods and distinguishes FAS eligible counts from sample sizes.
- [Complete-panel figure evidence bundle](results/PUBLICATION_COMPLETE_FIGURE_BUNDLE_20260926.md)
  verifies outside the checkout and includes reproducible retained-count arithmetic.
- [Eight-method OrthoBench strata](results/OB_COMPLETE_STRATA_RESULT_20260928.md)
  retain all 112 frozen method/bin rows and independent rational-arithmetic
  verification. The [three-metric figure](results/OB_COMPLETE_STRATA_FIGURE_20260928.md)
  retains all empty bins and adds no confidence intervals. This is separate
  from the SwissTrees bundle and the original three-method uncertainty panel.
  Its [direct-evidence bundle](results/OB_COMPLETE_FIGURE_BUNDLE_20260928.md)
  was verified outside the checkout with isolated standard-library Python;
  this does not provide standalone plotting or native inference.
- [Complete-panel relocated arithmetic check](results/SWISS_COMPLETE_PORTABLE_REPRODUCTION_20260926.md)
  reproduces all 24 endpoints from a three-file committed export using
  isolated Python; a file-access trace checks that original checkout paths
  were not accessed. This is not native inference or raw-source admission.
- [Complete descriptive strata](results/SWISS_COMPLETE_STRATA_20260926.md) and
  [eight-method feature figures](results/SWISS_COMPLETE_STRATA_FIGURES_20260926.md)
  retain the original bins and feature admissions; historical recipes below
  remain available with their original missing-method state.
- [Complete corrected eight-method table](results/qfo_corrected_comparison_20260926_v7/scores.md),
  [admitted recovered OrthoMCL scores](results/qfo_recovered_score_admission_22176.json),
  [complete SwissTrees uncertainty](results/qfo_recovered_swiss_uncertainty_22178.json),
  [24-endpoint arithmetic reproduction](results/qfo_recovered_swiss_reproduction_22178.json),
  and [figure/endpoint tables](results/corrected_swiss_comparison_figure_20260926/endpoints.md).
  The [job 22178 batch recipe](results/qfo_recovered_swiss_comparison_batch_20260926.sh)
  pins the recovered-score-compatible executor and comparison digest.
  The seven-method recipes below are retained historical stages, not the
  latest complete panel. All seven earlier contrasts remain unchanged.
- [Corrected figure supplements](results/PUBLICATION_CORRECTED_FIGURE_BUNDLE_20260919.md)
  include the seven-method FastOMA comparison exported at `27113f4`, with
  isolated relocated verification and retained-count arithmetic reproduction.
  Earlier factorial and composition-strata exports are preserved separately;
  none is a complete publication archive or a controlled timing admission.
- [Working manuscript](results/PUBLICATION_MANUSCRIPT_DRAFT_20260916.md)
  describes the frozen scientific configuration at `7f3a9e4`. Current
  development source is not interchangeable with that baseline.
- [Claim-to-evidence checklist](results/PUBLICATION_CLAIMS_20260916.md)
  distinguishes supported findings from limitations and missing evidence.
- [Corrected QfO protocol](results/QFO_CORRECTED_RELEASE_PROTOCOL_20260917.md)
  and [input validation](results/QFO_CORRECTED_INPUTS_STAGED_20260918.md)
  identify the corrected inputs. Historical QfO results remain separate.
- [Corrected seven-method table](results/qfo_corrected_comparison_20260923_v6/scores.md)
  contains admitted point estimates, not a complete eight-method ranking.
  The six-metric mean is a project-defined secondary summary, not a native
  QfO endpoint. FastOMA uses a supplied corrected OrthoFinder species tree;
  OrthoMCL corrected scores were unavailable in this historical snapshot,
  not in the current eight-method table above. The associated
  [paired SwissTrees results](results/qfo_fastoma_swiss_uncertainty_22098.json)
  retain all eight planned contrasts and all 24 endpoints for correction,
  with seven contrasts estimable. No general superiority is established.
- [OrthoBench paired analysis](results/ORTHOBENCH_UNCERTAINTY_20260916.md),
  [novel-taxa evaluation](results/YGOB_FROZEN_INTERPRETATION_20260916.md), and
  [Three Kingdoms audit](results/THREE_KINGDOMS_PAIR_COUNT_AUDIT_20260918.md)
  retain their distinct reference universes. Novel taxa do not imply
  family-disjoint validation; BUSCO recovery is not proteome-wide accuracy.

## Verified Search Diagnostics

- [Failed-query membership](results/qfo_recovered_failure_membership_20260926.json)
  verifies that all 53 failed OrthoMCL queries are absent from the recovered
  native graph index and final groups. [Reference exposure](results/qfo_recovered_reference_impact_20260926.json)
  measures direct overlap with benchmark references and annotations; it does
  not bound repaired-search effects on clustering or accuracy.
- [Residue-to-group trace](results/qfo_native_residue_group_trace_20260926.json)
  places the seven legacy residue-deletion proteins in three single-species
  groups with zero incident submitted cross-species pairs. It does not
  reconstruct a residue-preserving search.
- [Forced-candidate protocol](results/OB_FORCED_CANDIDATE_PROTOCOL_20260926.md)
  and [independent raw-output audit](results/ob_forced_candidates_audit_20260926.json)
  document 30,496 prefilter-excluded reference pairs passing unchanged HMM
  scoring. All 33,098 previously scored pairs retain exact scores and E-values.
  Reference-conditioned selection makes this a mechanism diagnostic, not
  unbiased recall, matched DIAMOND sensitivity or whole-pipeline F1 evidence.
- [Figure reproduction](results/FORCED_CANDIDATE_FIGURE_REPRODUCTION_20260926.md)
  provides the committed source export and isolated plotting command, verified
  outside the checkout in a fresh venv with a [hash-locked wheel set](results/forced_candidate_plot_linux_py310_requirements_20260926.txt).
  The PNG is byte-identical on the tested Linux x86_64/Python 3.10 platform.
  This plotting environment does not supply the native inference toolchain.

## Reproduce Statistics

The [complete eight-method OrthoBench replay](results/OB_COMPLETE_PORTABLE_STATISTICS_20260928.md)
works outside the checkout with one standalone checker, a 114 KB derived-count
file, Python and NumPy. It reproduces all 21 paired contrasts and seven family
win/tie/loss triples without historical absolute paths. Its receipt explicitly
distinguishes statistical replay from raw-input verification or inference.

The updated seven-method SwissTrees result was generated by job 22098 using
the unchanged executor at `10338e2a5046e522f2c1990e791532d3376e8982`.
Its [batch recipe](results/qfo_fastoma_swiss_comparison_batch_20260923.sh)
documents the exact source tree, inputs and protocol pins; invoking sbatch
would submit a job, unlike the arithmetic-only command below. Historical
six-method results remain retained. Later reviewed changes update FastOMA
admission paths and support recovered OrthoMCL admission/provenance without
changing the statistical kernel. The [new consolidation executor receipt](results/qfo_recovered_swiss_executor_20260923.json)
pins commit `0b8efa9c6b20b5bdec41f798c136b54132e7db37`, 123 passing tests and a
seven-method raw-count reconstruction with exactly matching estimates and
intervals. This executor is prepared for eventual recovered OrthoMCL evidence,
not evidence that corrected OrthoMCL scores exist. Follow the
[recovery handoff](results/QFO_RECOVERY_BPO_HANDOFF_20260923.md) only after actual
score admission. Do not substitute the new executor when reproducing the
historical executor's byte-level provenance.

Recheck the seven-method arithmetic from its retained audited family counts
in this workspace, where referenced provenance files are available:

```bash
python benchmark_tools/reproduce_corrected_swiss_comparison.py \
  --results benchmark_tools/results/qfo_fastoma_swiss_uncertainty_22098.json \
  --results-sha256 121cc8adf3cd63878f19006ec5500a13f879d042eccd565e1ae5f5a863f43fb1 \
  --output /tmp/orthohmm-fastoma-swiss-reproduction-new.json
```

The retained reproduction matches 21 available endpoints within `1e-12`;
the three OrthoMCL endpoints are missing, not zero. This rechecks arithmetic
and referenced hashes, not a fresh raw-count reconstruction or native run.
Repeating that full provenance check on another machine requires restoring
its referenced evidence paths. For a narrower portable arithmetic check,
the verifier at commit `5f3088a` can run with only the retained report and
NumPy; it does not follow those historical paths:

```bash
python -I -B benchmark_tools/reproduce_corrected_swiss_comparison.py \
  --retained-counts-only \
  --results benchmark_tools/results/qfo_fastoma_swiss_uncertainty_22098.json \
  --results-sha256 121cc8adf3cd63878f19006ec5500a13f879d042eccd565e1ae5f5a863f43fb1 \
  --output /tmp/orthohmm-swiss-portable-new.json
```

[Executed relocated evidence](results/swiss_portable_reproduction_20260923.json)
reproduces all 21 available endpoints within `1e-12` using Python 3.10.13
and NumPy 2.2.6. The isolated `/tmp` export contained the verifier, count
report, analysis requirements/lock and project license, all verified against
committed Git blobs. No sibling analysis modules or historical evidence
files were needed. This mode explicitly does **not** repeat raw-source
admission, inference or scoring; it is not the complete archival workflow.

Run from the repository root with Git history available, including the
explicit revisions below. Use a separate analysis interpreter matching
[analysis requirements](swiss_analysis_requirements.txt) and the associated
hash lock described in the linked execution records. The commands do not
create or install that environment. In this workspace its existing path is
`benchmarks/work/swiss_analysis_env_20260917/bin/python`; replace that path
on another machine. Use fresh output and report destinations on each run.

OrthoBench factorial, from 70 retained reference-family sufficient statistics:

```bash
python benchmark_tools/reproduce_orthobench_factorial.py \
  --revision a8a9f57734902715df7e9220cba2b3e0106fd5f3 \
  --python benchmarks/work/swiss_analysis_env_20260917/bin/python \
  --output /tmp/orthohmm-ob-factorial-reproduction-new \
  --report /tmp/orthohmm-ob-factorial-reproduction-new.json
```

[Executed evidence](results/ORTHOBENCH_FACTORIAL_REPRODUCTION_20260918.md):
eight cells, 20,000 paired draws, seed 20260918; all scientific fields
reproduced exactly outside the checkout.

Corrected QfO SwissTrees factorial, from 18 retained families:

```bash
python benchmark_tools/reproduce_corrected_factorial.py \
  --repo . --revision 0864f9b \
  --python benchmarks/work/swiss_analysis_env_20260917/bin/python \
  --output /tmp/corrected-factorial-reproduction-new \
  --report /tmp/corrected-factorial-reproduction-new.json
```

[Executed evidence](results/QFO_CORRECTED_FACTORIAL_REPRODUCTION_20260919.md):
eight cells, 100,000 paired draws, seed 20260922, 42 multiplicity-adjusted
endpoints; all numerical fields reproduced exactly. Additional executed
workflows cover the [sequence-search control](results/QFO_SEQUENCE_REPRODUCTION_20260918.md)
and [historical SwissTrees domain strata](results/SWISS_DOMAIN_RELOCATED_REPRODUCTION_20260917.md).
Do not relabel historical-release strata as corrected-release analyses.

These reproductions use committed counts, not historical absolute paths
embedded as provenance. They do not repeat native inference, conversion,
official scoring, raw-reference validation or independent biological testing.
The runners record installed dependencies but do not establish an OS-hermetic
environment or cross-platform equivalence.

## Native Execution And Timing

OrthoBench raw inputs and its upstream scorer are acquired separately:

```bash
git clone --no-checkout https://github.com/davidemms/Open_Orthobench.git /tmp/orthobench-source-new
git -C /tmp/orthobench-source-new checkout --detach 872d6f30592ab5ff837224db16a514b3f2bb916a
python benchmark_tools/verify_orthobench_acquisition.py \
  --checkout /tmp/orthobench-source-new --output /tmp/orthobench-acquisition-new.json
```

[Executed reacquisition evidence](results/orthobench_source_reacquisition_20260923.json)
verifies 95 files against the pinned upstream blobs. The 12 FASTAs, 70
reference groups and 11 low-certainty files also match the retained
scientific input manifests. The scorer and README are checked but not
executed. Raw files and upstream code remain outside the publication bundle;
these commands do not resolve redistribution rights or rerun inference.

YGOB inputs can be acquired directly from the provider without including raw
data in the publication bundle. From the repository root, use fresh paths:

```bash
python benchmark_tools/acquire_ygob_publication_inputs.py --output /tmp/ygob-source-new
python benchmark_tools/prepare_ygob_validation.py \
  --candidate /tmp/ygob-source-new \
  --audit benchmark_tools/results/ygob_overlap_20260916.json \
  --output /tmp/ygob-prepared-new --summary /tmp/ygob-prepared-new.json
```

The acquisition step uses only the standard library; preparation requires
the recorded Biopython environment. All three downloaded files must match
frozen sizes and SHA-256 values. Failed downloads are preserved and are not
silently retried. [Executed acquisition evidence](results/ygob_source_reacquisition_20260923.json)
matched the original inputs; regeneration matched all 16 prepared FASTAs
and the reference file byte-for-byte. HTTP transport is unauthenticated;
the hashes bind previously retained bytes. Raw YGOB redistribution remains
uncleared, and this does not repeat native inference or scoring.

Full native reproduction needs the exact input, software, runtime and
conversion records referenced by each experiment. The
[CPU baseline build](results/PUBLICATION_BASELINE_CPU_BUILD_20260919.md)
and [patched installer record](results/PUBLICATION_INSTALLER_SECURITY_20260919.md)
are bounded installation evidence, not complete phylogenetic reproduction.
Do not use the retained vulnerable historical installer lock for a new install.

The [DGX original-panel disposition](results/DGX_POSTRUN_ADMISSION_20260918.md)
admits native outputs but retains timings as descriptive only. The
[replacement specification](results/DGX_SCALING_REPLACEMENT_PROTOCOL_20260920.md)
and [long-run amendment](results/DGX_SCALING_LONG_RUN_AMENDMENT_20260920.md)
do not authorize a launch: environmental policy and executor integration
remain incomplete. Do not claim matched-resource speedups from the old panel.
The user deferred further DGX work on September 23. No further DGX access
requests or timing submissions are planned without renewed direction;
controlled comparative timing remains an unmet publication requirement.

## Archive And Outstanding Work

A [verified source-only baseline archive](results/PUBLICATION_FROZEN_SOURCE_ARCHIVE_20260920.md)
now preserves the 43 selected source files at scientific revision `7f3a9e4`,
with exact Git-blob/mode verification and repeat-export checks. It is a local
source component, not the analysis/data/runtime bundle or a public release.
Historical build behavior is preserved rather than silently replaced by
development packaging fixes.

A separate [setup-overlay installation](results/PUBLICATION_FROZEN_OVERLAY_INSTALL_20260926.md)
now builds and installs the frozen scientific sources with only the tested
packaging setup replaced. All 33 shipped scientific source files match frozen
Git blobs, and standard/high-sensitivity installed fixtures pass. This is a new
compiled artifact with its own lock, not the original benchmark runtime or
full phylogenetic reproduction; the source-only archive remains unchanged.

The [exact overlay artifact inventory](results/FROZEN_OVERLAY_ARTIFACT_INVENTORY_20260926.md)
binds all 11 installed wheel identities to provider notice metadata and verifies
all 42 RECORD entries, the embedded license and direct native linkage of the
new project wheel. Ten third-party wheels match the prior inventory; the new
OrthoHMM wheel is not silently treated as the earlier development artifact.
Full transitive dependencies and selected-archive redistribution review remain
separate requirements.

The subsequent [installed inferred-phylogeny fixture](results/PUBLICATION_FROZEN_PHYLOGENY_INSTALL_20260926.md)
executes HMM search, candidate expansion, marker-based species-tree inference,
gene-tree inference and reconciliation on 16 synthetic proteins. It verifies
fresh tree work, exact input coverage, pair integrity and parsed tree tips.
This closes the small-fixture pipeline check, not full publication-dataset
reproduction or packaging of external MAFFT/FastTree installations.

A [traced external-tool relocation](results/PUBLICATION_RELOCATED_PHYLOGENY_TOOLS_20260926.md)
copies MAFFT's entrypoint/helper files and FastTree locally, explicitly sets
`MAFFT_BINARIES`, and obtains byte-identical fixture outputs without traced
access to either original tool prefix. Python and system libraries are not
relocated, and binary redistribution/complete dependency closure remain open.

The [FastTree source acquisition](results/PUBLICATION_FASTTREE_ACQUISITION_20260926.md)
now pins six upstream release files, including source and notices, and matches
the installed executable byte-for-byte to upstream v2.2.0. Its executable
requires AVX2; this is not a portable source build or complete redistribution
review. The workflow downloads without execution and leaves installed tools
unchanged.

A separate [baseline-x86-64 FastTree source build](results/FASTTREE_BASELINE_BUILD_RESULT_20260926.md)
now produces byte-identical binaries in two private builds. An initial help-banner
validator failure is retained; after correcting that harness, the same binaries
pass the installed fixture and all four phylogeny readbacks, with five primary
outputs identical to the historical fixture. No old binary or benchmark runtime
was replaced. Full-dataset equivalence, execution on older CPUs and complete
dependency/redistribution review remain unproven.

A [MAFFT core source build](results/PUBLICATION_MAFFT_SOURCE_BUILD_20260926.md)
now reconstructs all 34 historical helper files byte-for-byte from the pinned
official 7.525 archive on this host. With its new prefix and explicit helper
path, the installed phylogeny fixture again produces identical groups, pairs
and selected trees. Optional RNA extensions were not built; source-archive
notices and full runtime/environment review remain separate requirements.

### Full Installed OrthoBench Readback

The [fresh installed run](results/INSTALLED_ORTHOBENCH_PROTOCOL_20260926.md)
was completed job 22179. It uses
the frozen scientific source with the packaging overlay and rebuilt MAFFT.
It is shared-host reproduction, not controlled comparative timing.

Independent readers have passed on the installed synthetic fixture and retained
historical OrthoBench p1_c1_r1: [trees and pair structure](results/PHYLOGENY_STRUCTURE_READBACK_20260926.md),
[sequence content and supermatrix](results/PHYLOGENY_SEQUENCE_READBACK_20260926.md),
[events, membership filtering and complete pair sets](results/PHYLOGENY_EVENT_READBACK_20260926.md),
and [hierarchy/selection](results/PHYLOGENY_HIERARCHY_READBACK_20260926.md).
The hierarchy is pre-membership-filtering; benchmarked final root groups and
native pairs are post-filtering. None is silently substituted for another.

Job 22179 subsequently completed successfully. All five commands below passed;
the [readback](results/INSTALLED_ORTHOBENCH_READBACK_20260926.md) records a
0.284504-percentage-point F1 decrease and a different partition. Historical
scores remain unchanged; this particular installed run did not reproduce them.
The subsequent controlled recovery is documented in the next section.

For a successful native completion, run the following from the repository root
in the audit environment (tested with Python 3.10.13, DendroPy 5.0.8 and
Biopython 1.86). Existing receipts for this run must not be overwritten;
these commands document their creation. Each reader refuses to overwrite a receipt. Stop on any
failure and preserve it; do not automatically resume or restart inference.

```bash
set -e
run=benchmarks/work/publication_installed_orthobench_20260926
phylo="$run/inference/orthohmm_phylogeny"
python -m benchmark_tools.audit_installed_orthobench \
  --repo . --directory "$run" --job 22179 --output "$run/readback_scores.json"
python -m benchmark_tools.audit_phylogeny_structure \
  --directory "$phylo" --input "$run/input" --output "$run/readback_structure.json"
python -m benchmark_tools.audit_phylogeny_sequences \
  --directory "$phylo" --structure "$run/readback_structure.json" \
  --output "$run/readback_sequences.json"
python -m benchmark_tools.audit_phylogeny_events \
  --directory "$phylo" --structure "$run/readback_structure.json" \
  --constraints "$run/inference/orthohmm_working_res/phylogeny_candidate_merges.json" \
  --output "$run/readback_events.json"
python -m benchmark_tools.audit_phylogeny_hierarchy \
  --directory "$phylo" --events "$run/readback_events.json" \
  --output "$run/readback_hierarchy.json"
```

The first command gates native/scheduler success and frozen input/source/tool
identity, scores all reference families and compares partitions without group
label dependence. Explicitly review every score/partition difference and all
receipts before any reproduction claim. The companion readers establish
output/rule consistency conditional on saved trees; they do not prove optimal
rooting, biological truth or correctness of the HMM search. They do not by
themselves establish publication readiness or overwrite historical scores.

### Tested Canonical OrthoBench Recovery

The [completed recovery](results/CANONICAL_OB_PHYLOGENY_RESULT_22324.md)
reproduces all 59,770 historical root groups and all 70 reference-family score
records, including F1 74.10607351873405%. This is a same-host, stage-composed
reproduction, not a claim that the clean wheel alone reproduces historical
outputs or that a portable one-command publication environment is complete.

| Stage | Retained execution | Evidence and scope |
| --- | --- | --- |
| Frozen scientific sources with setup-only packaging overlay | Full installed job 22179 | [Installed readback](results/INSTALLED_ORTHOBENCH_READBACK_20260926.md); scientific revision `7f3a9e4`, independent input/output audits; preserve its differing F1 73.821569% |
| Fresh search evidence | Checkpoint from 22179 | [Search comparison](results/INSTALLED_OB_SEARCH_COMPARISON_20260926.md); same directed nonself hit keys with tiny score differences; do not normalize the checkpoint again |
| Clustering/profile replay | Job 22320, private Leiden 0.11 arm | [Downstream readback](results/OB_DEPENDENCY_REPLAY_RESULT_22320.md); identical historical pre-candidate partitions, 62,885 refined seed groups |
| Canonical candidate formation | Job 22323, `fresh_full_self_control` | [Five-arm experiment](results/OB_CANONICAL_CANDIDATE_RESULT_22323.md); sort by directed query/target indices without changing scores, thresholds or self-hit policy; 54,445 historical candidate families |
| Membership constraint binding | Independent semantic audit | [Constraint readback](results/CANONICAL_MEMBERSHIP_CONSTRAINTS_20260926.md); 8,440 directed constraints identical to historical consumer inputs and order |
| Fresh species/gene trees and reconciliation | Job 22324 | [Frozen protocol](results/CANONICAL_OB_PHYLOGENY_PROTOCOL_20260926.md); 32 CPUs, inferred species tree, minimum-variance rooting, species-overlap root rule, positive-paralogy pairs, membership constraints, no old tree checkpoints |
| Native admission, four scientific readers and frozen scoring | Job 22325 | [Workflow](results/CANONICAL_OB_PHYLOGENY_ADMISSION_20260926.md); exact partition and full score-object agreement after successful terminal accounting |

The private runtime used the retained Leiden 0.11 distribution, including its
bundled native libraries, with NumPy 2.2.6 and Python igraph 1.0.0. A version
string alone is not an adequate replacement for the distribution/file pins.
The clean overlay wheel lock includes Leiden 0.12; do not silently substitute
that lock for the recovery environment. See the
[private distribution experiment](results/OB_LEIDEN_OVERLAY_PROBE_20260926.md).
The [verified recovery wheel](results/LEIDEN_RECOVERY_WHEEL_20260926.md)
now supplies a provider-digest-checked artifact whose 15 payload files match
that private distribution exactly. Its single-component hash lock is not a
complete environment. A [fresh 11-wheel recovery installation](results/PUBLICATION_RECOVERY_INSTALL_20260926.md)
now passes source/package byte checks and standard/high-sensitivity fixtures,
using a separate complete Python lock. Full native phylogeny and benchmark
reproduction in that new installation have not yet been validated.
The canonical-order policy remains benchmark-only, not an enabled production
default. Do not present this as a new accuracy-tuned method or transfer the
recovery result to QfO or another platform without validation.

Native phylogeny used rebuilt MAFFT 7.525 and the retained AVX2 FastTree 2.2
binary. The separate baseline-x86-64 FastTree build was tested on a fixture,
not substituted into this full-data recovery. Exact interpreter, environment,
source, tool and input identities are in the local plan
`benchmarks/work/canonical_ob_phylogeny_20260926/plan.json`, SHA256
`2ded17c15cf3510e02f1683973edf211dc8b7f7f86037dcf83e121a948dfafdc`.
The [native submission](results/canonical_ob_phylogeny_submission_22324.json)
and [audit submission](results/canonical_ob_readback_submission_22325.json)
retain the executable commands and source/runtime preflight. They are evidence
of past executions, not portable launchers that acquire all missing artifacts.

To re-audit this completed local experiment without launching inference, choose
an output directory that does not exist:

```bash
python -m benchmark_tools.readback_canonical_ob_phylogeny \
  --repo . \
  --directory benchmarks/work/canonical_ob_phylogeny_20260926 \
  --job 22324 \
  --output /absolute/fresh/canonical_ob_readback
```

This requires the retained plans, absolute-path inputs, native outputs,
historical references, scheduler accounting and tested audit dependencies. It
is not a relocated native reproduction. The original `readback` directory must
not be overwritten. All audit stages must succeed before claiming a reproduced
score; a partial receipt or process exit alone is insufficient.

Fresh search/candidate artifacts were reused, so the 1,555.40-second replay is
phylogeny-stage timing only. Shared-host resources are descriptive. Remaining
integration work includes packaging the exact dependency distribution, exposing
a tested full-pipeline reproduction entrypoint with explicit ordering policy,
and validating its data/tool restoration without relying on workstation paths.
Do not change the frozen scientific baseline or silently promote the ordering
policy while doing that integration.

The [updated OrthoBench provenance register](results/OB_PROVENANCE_REGISTER_20260927.md)
and [QfO OrthoFinder consolidation](results/QFO_ORTHOFINDER_PROVENANCE_CONSOLIDATED_20260926.md)
distinguish available historical inference, conversion and scoring evidence.
The QfO sequence-only row is a full-run checkpoint diagnostic, not a separately
timed sequence-only inference. The OrthoHMM high-sensitivity cached replay
is 319.467197 seconds, not a full inference measurement; the historical
phylogenetic full run is 3274.102675 seconds, not the fresh installed recovery.
These scopes and memory accounting differ. Other cross-dataset provenance
remains incomplete.

The [rights register](results/PUBLICATION_DATA_RIGHTS_20260918.md) separates
project code, reference datasets, upstream software and binaries. Local
availability or inclusion in an evidence bundle is not blanket permission
to redistribute. Original TreeFam trees and mapping remain
[unrecovered](results/TREEFAM_SOURCE_RETRIEVAL_20260918.md); the pooled reference
does not supply original family identities for uncertainty estimation.

### Remaining Publication Gates

- Controlled timing: the 27 replacement scaling runs remain unexecuted.
  The local Threadripper `bizon` is the approved host; DGX access or permissions
  are not prerequisites. Complete the accounting/environmental validation and
  establish a verified quiet window before launch. Do not substitute
  uncontrolled shared-host timings or restart DGX work.
- Uncertainty: native VGNC, GO/EC, FAS and the secondary summary lack admitted
  paired intervals. The [rare-event screen](results/SPARSE_DYADIC_F1_RESULT_20260926.md)
  rejects blanket promotion of the tested Wald approach. More implementation
  alone cannot justify its sampling/dependence assumptions. Original TreeFam
  source families/mapping remain an external acquisition gap; no inferred
  connected-component family labels may be substituted.
- Generalization: the frozen YGOB transfer supports a bounded novel-taxon
  evaluation, not family-disjoint or universal generalization. No outcome-driven
  method changes are authorized by that evidence; such changes require a new
  independent confirmation.
- Executable packaging: inference, reader and statistical components have
  bounded relocation checks, but their combined acquisition/install/run/score
  workflow and full-data restoration are not yet a complete portable release.
  A final bundle must specify which historical source versions each command
  consumes and preserve both failures and successful executions.
- Distribution: OrthoBench/YGOB remain acquisition-only. Final selected-file
  notices, binary/container dependencies, dataset exclusions and source
  obligations need review; local archives are not redistribution clearance.

These are open requirements, not claims that the study is ready except for
typesetting. No public archive DOI, external
submission or completed publication release is asserted. The
[progress ledger](results/PUBLICATION_PROGRESS.md) records dated actions;
use live scheduler state, not historical ledger prose, to decide whether
a job is still running. Preserve failed and incomplete experiments rather
than silently rerunning or selecting favorable outcomes.
