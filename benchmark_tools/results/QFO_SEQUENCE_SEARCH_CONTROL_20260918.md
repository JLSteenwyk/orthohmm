# Corrected QfO Sequence-Search Control

## Scope

Extend the exploratory OrthoBench initial-search replacement to the frozen
corrected QfO2020_04 panel:984,137 proteins in78 species. This remains
development-exposed evidence, not independent confirmation. No defaults
change, and no QfO endpoint is used to select parameters.

Use the existing DIAMOND2.1.11 binary and unchanged search_command from
prepare_sequence_search_control.py:32 threads, very-sensitive, E-value1e-4,
BLOSUM62, gap11/1, composition1, masking1, all targets and one HSP. Query
all proteins against each species-specific database. Retain direction and
self hits; apply no identity/coverage cut. Preserve raw scores, bit scores,
E-values and lengths for audited normalization once by sqrt(query*target
length). Original OrthoBench manifest supplies binary/settings provenance
only; no original hit table or checkpoint is reused.

Primary comparison is P0/C0/R0 downstream graph processing versus the
corrected initial-HMM-search control, with the same frozen graph settings.
Retain all-hit and post-search top100 diagnostic variants as specified in
SEQUENCE_SEARCH_CONTROL_PROTOCOL_20260916.md. The post-search cap is not
the HMM prefilter cap. Equal E-values and identical graph parameters do
not establish equal sensitivity, score calibration or computational cost.
Do not claim that this alone solves the matched-sensitivity requirement.

## Prepared Inputs

Preparation implementation3cbe3a3 verifies pinned corrected staging and
direct inventory, every FASTA hash, exact gene/species scope and the existing
DIAMOND binary hash/version. It writes no search result and authorizes no
execution.24 focused preparation/conversion tests passed in0.26s.

Root: benchmarks/work/qfo_sequence_search_control_v1.

| Artifact | Bytes | SHA-256 |
| --- | ---: | --- |
| manifest.json | 205271 | fcef19d9c745c5806197aa219dec1c89bd5d90082a35f24dadafa3422f461436 |
| queries.fasta | 584039577 | 00e106d9e6f0eaf44cbaca0fbde6c63e2f587b848e6ca9cbaeb5d7a2268b5519 |
| gene_metadata.json | 80596465 | 3ce892e98e5254270729d44c187ac57fb706154385d590e3048a98c28f49b511 |

Post-preparation checks rehashed source, queries, metadata and binary;
all78 target directories are empty. Large derived files remain untracked.

## Remaining Gates

The [graph memory planning review](QFO_GRAPH_MEMORY_PLANNING_20260918.md)
documents full-length allocations in the frozen reciprocal-hit implementation
and a tested, checkpoint-bound array-payload diagnostic. Its completed
OrthoBench calculations are not QfO peak-RAM predictions. Apply it only to
admitted corrected checkpoints before selecting graph-replay resources;
downstream graph/clustering costs remain outside its one-snapshot bounds.

- Freeze a corrected-only runner with pre/post checks and full-panel
  success accounting; retain failure logs and forbid partial inference.
- Measure database and search phases separately. Shared-workstation cost
  is descriptive, not dedicated timing evidence. Do not load the active DGX.
- Validate hit ID ownership, exact lengths, unique ordered pairs, finite
  scores and target scope before numeric checkpoints. SQLite/memmap keeps
  conversion chunked but does not prove the downstream graph fits memory.
- Review observed output size and peak memory before graph replay; do not
  truncate all-hit data silently to meet resource limits. The historical
  OrthoBench all-hit checkpoint has100,099,147 hits, not a prediction of
  corrected QfO output size.
- Compare hit overlap/coverage and score distributions against the admitted
  corrected HMM checkpoint before accuracy interpretation, without fitting
  scale factors to benchmark outcomes. Independently validate and score
  both complete diagnostic variants, retaining negative findings.

## Runner Submission

Runner b66d3f374224b89f782c98a4b29a9c33b6cefdc7 is frozen at
benchmarks/work/publication_qfo_sequence_search_v1.37 focused tests pass;
actual full-input preflight succeeds without creating execution output.
Job21789 was submitted2026-09-18T10:57:20 with32CPUs/192GiB/7days onbizon,
no requeue. Scheduler confirms the requested allocation and pending state.
The seven-day limit is an execution allowance, not a runtime prediction.

All78 targets run sequentially, measuring makedb/search separately and
retaining failed phase evidence. Full postflight rechecks input/source/tool
identities and every hit/database/log/timing record before marking the
panel complete_pending_numeric_validation. That status is not accuracy
admission. The preparation manifest remains unchanged; its false execution
authorization is a preparation-state field, not a completed-run claim.
Numeric conversion, graph replay, hit diagnostics and scoring remain open.

Independent panel admission21790 follows afterany:21789, using frozen
6cf2e92074733cab036c1dd1674467bd37ce2313 at
publication_qfo_sequence_search_admission_v1. It requires successful
terminal accounting and verifies all78 target commands, phase outcomes and
retained file identities.2CPU/64GiB/24h/bizon/no-requeue; submitted
2026-09-18T11:01:31.55 focused tests and actual live-job rejection pass.
Future output: benchmarks/work/qfo_sequence_search_admission_20260918.json.
This is execution admission only; numeric hit validation is still required.

Numeric conversion21791 follows afterok:21790 with frozen
f5c23b81e3cf72573fafd5d1bb087fd2e61392d2 at
publication_qfo_sequence_numeric_v1. Submitted2026-09-18T11:05:46,
2CPU/192GiB/7days/bizon/no-requeue. It checks metadata directly against
corrected FASTAs, ingests all78 outputs with duplicate-pair rejection and
writes both all-hit and post-search top100 audited numeric checkpoints to
benchmarks/results/qfo_sequence_numeric_v1.46 focused tests pass. No
production checkpoint is available yet, and frozen graph replay plus
independent checkpoint equivalence/scoring checks remain required.

Independent numeric admission21792 follows afterany:21791, frozen at
617588f0f0acef0134b3e5fe18480e00eaa22129 in
publication_qfo_sequence_numeric_admission_v1. Submitted2026-09-18T11:21:33
with2CPU/192GiB/7days/bizon/no-requeue. It independently reconstructs source
hits in a separate SQLite database and compares exact normalized tuples,
gene/species maps, and streaming per-query/species top100 ranks against both
hashed checkpoints. Terminal accounting, frozen converter/admitter sources,
input identities and pre/post file hashes are required.64 focused tests
pass, and actual pending-conversion preflight refuses without output.
Future report: benchmarks/work/qfo_sequence_numeric_admission_20260918/report.json.
No production equivalence result exists yet. This gate does not authorize a
claim of matched search sensitivity, scientific accuracy or graph feasibility.

Label-free comparison21793 follows afterok:21720:21792, requiring both the
native corrected HMM admission and independent DIAMOND numeric admission.
Executor3a83f9db2981c7bad3d984738c599cf0d1bb4585 is frozen at
publication_qfo_hit_coverage_v1. Submitted2026-09-18T11:33:34 with
2CPU/192GiB/7days/bizon/no-requeue. A separate disk-backed HMM copy sorts
native species-pair batches without changing gene IDs or score values;
duplicates fail rather than being collapsed. All-hit/top100 remain unchanged.
The comparison reports exact directed overlap and coverage plus fixed score
histograms; it does not estimate biological sensitivity or match calibration.
Future report: benchmarks/work/qfo_corrected_hit_coverage_20260918/report.json.
83 focused tests pass; actual pending-admission preflight refuses without
output. No production diagnostic or graph-feasibility result exists yet.

## SwissTrees Source-Bound Count Audit

`audit_qfo_sequence_swiss.py` accepts `--inventory`, `--inventory-sha256`,
`--baseline` and a fresh `--output`. The inventory contains a `variants` list
in exact order `p0_c0_r0`, `all_hits`, `top100`; each entry contains `variant`
and an `admission` record with absolute `path`, `bytes` and `sha256`.
The baseline is the pinned original SwissTrees count audit with SHA-256
`546bb5bd6957c8ea990324b79fa31f22b0ed721bc7d6b94b609ab15258f97183`.
Only its reference identities/labels/membership are reused, never predictions.

The auditor checks corrected admission identities, completed scorer/converter
accounting, embedded conversion equality, inventoried raw output and file
hashes before and after reconstruction. It compares exact pair labels and
members to the audited reference, verifies all native family scores and macro
statistics, and requires18 families/10765 relations. Output status is
`corrected_qfo_sequence_swiss_counts_verified`; uncertainty remains unadmitted.
It relies on upstream inference/scoring admissions rather than rerunning those
entire workflows.74 focused tests pass, including23 new audit tests and an
integrated synthetic file-bound audit with a post-admission tampering check.

No production count audit has run: all three corrected assessments are still
required. The numerical protocol is frozen in
`QFO_SEQUENCE_UNCERTAINTY_PROTOCOL_20260918.md`; its production bootstrap
runner, binding exact protocol/auditor/helper sources, remains to be added.

### Production Uncertainty Runner Added

`run_qfo_sequence_uncertainty.py` now accepts `--counts`, `--counts-sha256`,
`--protocol`, `--output` and `--markdown`. Supply the count auditor's JSON and
its independently recorded SHA-256, the unchanged uncertainty protocol above,
and two distinct fresh output paths. Run from the frozen source checkout used
for count assembly: source record paths must match as well as content hashes.

The runner pins the protocol and twelve source/helper files, checks all count
evidence records, calls the source-bound auditor again, and requires exact
equality with the saved counts before resampling. All recorded inputs are
rechecked afterward. Production controls are fixed at100000 replicates,
PCG64 seed20260923 and six adjusted endpoints; CLI overrides are not exposed.
Successful output admits only these conditional SwissTrees intervals and
keeps publication_ready=False. The runner does not rerun upstream inference
or native scoring and is not a hermetic software-environment installer.

55 focused tests pass, including a synthetic100000-draw invocation and
rejection of changed protocols, helpers, counts, evidence and output collisions.
No real corrected sequence-control intervals have been produced. The frozen
protocol's outstanding-work paragraph records its original freeze-time state;
this operational update does not alter its endpoints or statistical controls.
