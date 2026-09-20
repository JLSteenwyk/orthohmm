# Publication Reproduction Guide

Status: 20 September 2026, incomplete working package. This guide routes
readers to executed workflows and their evidence; it is not an archival
release or an assertion that every requirement is complete. Commands below
are statistical reproduction only and do not submit scheduler jobs.

## Method And Results

- [Working manuscript](results/PUBLICATION_MANUSCRIPT_DRAFT_20260916.md)
  describes the frozen scientific configuration at `7f3a9e4`. Current
  development source is not interchangeable with that baseline.
- [Claim-to-evidence checklist](results/PUBLICATION_CLAIMS_20260916.md)
  distinguishes supported findings from limitations and missing evidence.
- [Corrected QfO protocol](results/QFO_CORRECTED_RELEASE_PROTOCOL_20260917.md)
  and [input validation](results/QFO_CORRECTED_INPUTS_STAGED_20260918.md)
  identify the corrected inputs. Historical QfO results remain separate.
- [Corrected six-method table](results/qfo_corrected_comparison_20260919_v5/scores.md)
  contains admitted point estimates, not a complete eight-method ranking.
  The six-metric mean is a project-defined secondary summary, not a native
  QfO endpoint. FastOMA and OrthoMCL corrected scores are still pending.
- [OrthoBench paired analysis](results/ORTHOBENCH_UNCERTAINTY_20260916.md),
  [novel-taxa evaluation](results/YGOB_FROZEN_INTERPRETATION_20260916.md), and
  [Three Kingdoms audit](results/THREE_KINGDOMS_PAIR_COUNT_AUDIT_20260918.md)
  retain their distinct reference universes. Novel taxa do not imply
  family-disjoint validation; BUSCO recovery is not proteome-wide accuracy.

## Reproduce Statistics

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

## Archive And Outstanding Work

The [rights register](results/PUBLICATION_DATA_RIGHTS_20260918.md) separates
project code, reference datasets, upstream software and binaries. Local
availability or inclusion in an evidence bundle is not blanket permission
to redistribute. Original TreeFam trees and mapping remain
[unrecovered](results/TREEFAM_SOURCE_RETRIEVAL_20260918.md); the pooled reference
does not supply original family identities for uncertainty estimation.

Pending comparator and robustness analyses, controlled resource evidence,
remaining uncertainty limitations, final manuscript preparation and
file-level archive review remain open. No public archive DOI, external
submission or completed publication release is asserted. The
[progress ledger](results/PUBLICATION_PROGRESS.md) records dated actions;
use live scheduler state, not historical ledger prose, to decide whether
a job is still running. Preserve failed and incomplete experiments rather
than silently rerunning or selecting favorable outcomes.
