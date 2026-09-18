# Audited OrthoMCL BPO Checkpoint

`prepare_orthomcl_bpo_checkpoint.py` now connects the previously tested
streaming converter, complete source-to-BPO content auditor, native index
construction and independent native index validator. This is a preparation
component, not permission to execute inference from an unadmitted search.

The workflow requires a fresh directory, verifies the retained native runtime
and system-helper inventories, records input/helper hashes, converts with
the fixed legacy cutoff 1e-5, and checks every retained record against the
source HSPs. The guarded Perl launcher then invokes native OrthoMCL index
construction and checks every byte offset, EOF sentinel and inclusive query
range. Source/helper hashes and runtime inventories are rechecked before a
successful report is written. Native stderr is retained and requires review
even with exit code zero. Partial files and a failure report are preserved;
the directory cannot be implicitly resumed or overwritten.

## Retained Fixture

The actual installed-runtime run at
`benchmarks/work/orthomcl_bpo_checkpoint_probe_20260918/` used the previous
guarded parity fixture: four input proteins, eight HSP rows, seven contiguous
pair blocks, one cutoff-excluded block, six BPO records, three query ranges
and seven offset entries including EOF. Its report status is
`bpo_checkpoint_content_and_indexes_verified`.

Retained report: `orthomcl_bpo_checkpoint_probe_20260918.json`, SHA-256
`09f3d847efbb5f5e896f8b518cc9a458eaaf3092507b9e7ee0382655ad0fe119`.
This fixture exercises the installed native constructor, not a replacement
index implementation. The orchestration unit tests use explicit mocks;
the opt-in installed-runtime test and retained run provide additional evidence.

61 focused tests passed with `ORTHOHMM_LEGACY_BLAST_SMOKE=1`, including source
content failures, stage failures, diagnostics, altered source bytes, incorrect
index summaries, existing-directory and dangling-symlink rejection, native
BPO parity, and damaged native index rejection.
The full unit suite subsequently passed: 3,726 tests in 71.05 seconds with
the legacy-runtime tests enabled. All 18 retained input/helper/output file
records were rechecked successfully after the fixture run.

## Boundaries

- Search admission remains false in the checkpoint report. Source database
  identity, scheduler completion and sequence-specific BLAST diagnostics
  must be verified by the corrected search admission workflow.
- No corrected production BPO or inference was run. The future execution
  wrapper must bind admitted search records, Python runtime, scheduler
  resources and this source version before starting full preparation.
- Index construction and validation keep native offset/range data in memory;
  this is not constant-memory indexing or a timing benchmark.
- Conversion/index correctness does not establish final-group accuracy.
  Native inference, final-group admission and QfO scoring remain required.
