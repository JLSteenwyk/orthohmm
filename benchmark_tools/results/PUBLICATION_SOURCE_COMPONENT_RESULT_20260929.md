# Scientific And Workflow Source Component

The [executed packaging receipt](publication_source_components_20260929.json)
records the first combined source-only component, retaining scientific and
workflow sources in separate directories rather than upgrading the baseline.

| Component | Revision | Files |
| --- | --- | ---: |
| Scientific package | `7f3a9e40dd7e79f842cc2c11fb8b548f9a802806` | 43 |
| Benchmark workflow and unit-test source, licenses and two guides | `88a379345988f84d1b54b404fadc0c06eca7f418` | 1,664 |

The 1,707 payload files total 9,207,865 bytes; all 1,693 Python files pass
syntax compilation without execution or imports. Scientific source is the
historical 0.5.0 package; the workflow revision is deliberately different.
The exporter reads committed Git blobs, so unrelated dirty sample files and
uncommitted source do not enter this component.

Local archive:
`benchmarks/work/publication_source_components_20260929/orthohmm-0.5.0-scientific-7f3a9e4-workflows-88a37934.tar.gz`

- Archive size: 2,341,862 bytes.
- Archive SHA-256: `683067d6eb0cdb1ebd79ee089992556139b3b0e7c6767c6070c044ead38d1929`.
- `SOURCE_INDEX.json` size: 689,511 bytes.
- Index SHA-256: `f15c343cfb52442d77d1c58ba7292addbe3f190377b7ffe9bdbbd555b122d0b8`.

The archive was extracted to a fresh temporary directory and its bundled
verifier executed under `/usr/bin/python3 -I -B`, with a PATH containing no
Git executable. The complete result matches the original verification.
The temporary copy was removed afterward; the original component and archive
remain retained. The verifier neither imports project dependencies nor reads
the original source checkout or historical workstation paths.

## Reproduction

From a clone containing both revisions, select a fresh output:

```bash
python -B benchmark_tools/bundle_publication_source.py build \
  --repo . --revision 88a379345988f84d1b54b404fadc0c06eca7f418 \
  --output /absolute/fresh/publication-source-component
```

After transferring/extracting the retained archive, verify the external index
digest as well as the archive digest before trusting its contents:

```bash
python3 -I -B /relocated/component/workflow/benchmark_tools/bundle_publication_source.py \
  verify /relocated/component \
  --manifest-sha256 f15c343cfb52442d77d1c58ba7292addbe3f190377b7ffe9bdbbd555b122d0b8
```

35 focused source-component/frozen-archive tests pass. They include immutable
Git export from a dirty checkout, distinct source revisions, payload/mode
corruption, missing/extra/symlink paths, invalid mappings and scope flags,
coordinated file/index changes rejected by the external digest, existing-output
refusal, and isolated verification after deleting the synthetic source checkout.
The 22 warnings are the already-retained invalid escape strings in frozen
parser/writer source; those scientific files were not changed.

## Publication Boundary

This closes the previous source-component omission of benchmark workflow and
unit-test code. It does not close the full executable-study or public-release
requirements. Raw datasets/references, test sample datasets, predictions,
scoring outputs, result receipts/plans, figures/manuscripts, third-party source
archives, wheels, binaries and runtime/OS images are excluded. Such materials
need the separately documented acquisition/evidence/runtime components.

The included reproduction guide still links to excluded evidence files, and
test source alone does not ensure all tests run without fixtures/dependencies.
Historical scripts retain DGX commands and absolute paths for provenance, not
authorization. The local Threadripper remains the approved timing host.
No new benchmark, numerical-equivalence check, compiled-library attribution,
complete redistribution clearance, external deposition, Git tag or package
version is asserted. The archive is local, not a publication-readiness claim.

All 822 existing source pins for pending diagnostic job 22378 were independently
rechecked and still match. The exporter is a new module, not a modification of
that frozen set. No inference was restarted and no service/scheduler policy
or unrelated work was changed.
