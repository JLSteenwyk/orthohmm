# SwissTrees Raw Source Relocation

The raw duplication and historical-fragment exporters require retained files
outside main Git, not just derived tables. Provide an explicit, opt-in staging
and relocation workflow rather than bypassing these dependencies. The default
exporter mode still uses original paths and all existing source/statistic gates.
No scientific settings, raw manifests, counts, score tables or old archives change.

## Workflow

`relocate_swiss_raw_sources` checks a frozen source manifest and every raw
record before copying. Stage each unique content digest once, with one-MiB
streaming copies, check the copies and original sources again, and preserve
every original record occurrence in an immutable relative-path binding file.
Destination directories must be new and without symlink parents. Failed partial
attempts remain available for diagnosis; do not silently reuse or overwrite them.

Both exporters accept optional `--source-bindings` and
`--source-bindings-sha256` together. Require an independently retained binding
digest, the exact frozen manifest content identity, and the complete ordered
source list including duplicates. No missing/extra/reordered/repinned records,
arbitrary paths, symlink escapes or rights claims are accepted. Validate actual
restored bytes before computing rows and again before writing the success
manifest. New output provenance records the binding/helper identities and
the unchanged historical source records. Original manifests are not rewritten.

On the original source host, from the repository root:

```bash
python -m benchmark_tools.relocate_swiss_raw_sources \
  --source benchmark_tools/results/swiss_duplication_features_v2_20260923.json \
  --source-sha256 97b0c4755d6a9df258d5c3f60fc0d5d25f1e5c09c42216c754a245a67d1942ec \
  --records-key checked_inputs --output /absolute/new/duplication-inputs

python -m benchmark_tools.relocate_swiss_raw_sources \
  --source benchmark_tools/results/swiss_historical_fragment_admission_22117.json \
  --source-sha256 a480f32666bb96de3395213c7cad024c170318ca110c389f15e4756f8052adf8 \
  --records-key records --output /absolute/new/fragment-inputs
```

Retain the printed binding digest independently before transferring the complete
directory. Do not derive an expected digest from an untrusted received manifest.
On a restored checkout with the required analysis dependencies:

```bash
python -m benchmark_tools.export_swiss_duplication_strata \
  --counts benchmark_tools/results/qfo_recovered_swiss_uncertainty_22178.json \
  --counts-sha256 a599749d66433ec211ed1ab0e3a6a2a75eb99cb0dc57abc47c2d267930cbbe34 \
  --features benchmark_tools/results/swiss_duplication_features_v2_20260923.json \
  --source-bindings /restored/duplication-inputs/bindings.json \
  --source-bindings-sha256 2009be32157feddd2d0f6f4191bcbf5e0d33416b24e9a50e402f6aa023df8e58 \
  --output /absolute/new/duplication-table

python -m benchmark_tools.export_swiss_fragment_strata \
  --counts benchmark_tools/results/qfo_recovered_swiss_uncertainty_22178.json \
  --counts-sha256 a599749d66433ec211ed1ab0e3a6a2a75eb99cb0dc57abc47c2d267930cbbe34 \
  --admission benchmark_tools/results/swiss_historical_fragment_admission_22117.json \
  --admission-sha a480f32666bb96de3395213c7cad024c170318ca110c389f15e4756f8052adf8 \
  --source-bindings /restored/fragment-inputs/bindings.json \
  --source-bindings-sha256 0319644fbb580d87ab123676a7220f89f053cac89e75a5674c4295dd4d3f72f1 \
  --output /absolute/new/fragment-table
```

Those binding digests identify the retained local staging directories under
`benchmarks/work/swiss_raw_sources_{duplication,fragment}_20261002`, not a
public downloadable distribution. They are valid only for these exact manifests;
a newly staged binding from a differently located source manifest has its own
digest even if all raw content digests match. The relative artifact paths make
an unchanged staged directory portable after transfer.

## Verification

Initial 40 relocation cases pass in 0.40s; the overlapping first exporter
panel has 62 passes in 4.11s. Add destination-symlink and during-copy/validation
mutation checks. **187 focused cases pass in 8.67s**, zero failures/errors/skips,
including 43 relocation cases and all prior raw-exporter/count-selection tests.

The actual copied-tree test stages 1,823 project/data files (12,231,453 bytes)
and both raw directories in a new temporary location. An isolated Python
3.10.13 child executes both production exporters with explicit bindings:
all 32 duplication and 56 fragment rows equal their retained manifests;
both TSV and Markdown files are byte-for-byte identical to the 20260926 tables.
All 18/1,781 input records are checked, including 11/1,774 raw occurrences.
The fragment list contains one duplicated historical path; all 1,774 occurrences
remain, with 1,773 content files/606,744,832 bytes. Duplication has 11 content
files/93,699,174 bytes. Twenty-two imported project origins are staged, the
original-checkout canary is blocked, later original Python open events are zero,
and child subprocesses are forbidden. Temporary copied staging is removed;
the private source directories are retained locally. The last three unit tests
were added after that copied export; production code used there is unchanged.

[Machine-readable source and receipt pins](swiss_raw_source_relocation_20261002.json)
retain this evidence. Python-event guarding is not OS containment. Hash-checking
the native image and annotation evidence does not execute Darwin, recollect
UniProt annotations, rederive duplication features, rescore raw predictions,
regenerate bootstrap intervals or reproduce inference. This is exact raw-file
identity-gated table export, not a new biological admission or native audit.

## Remaining Boundaries

Preceding source-ee2a7af0 CI run 36992558116 was observed live. Docs, Linux
diagnostics and wheel succeed; macOS 3.10/3.11 fail while full/3.12/3.13 remain
live at that observation. This new relocation workflow is not in that source.
Do not restart jobs or infer sibling counts/new-patch confirmation. Existing
CI tests still require their raw inputs in default mode: this opt-in workflow
does not automatically provision them or claim full CI success.

At the later terminal observation, all five macOS test jobs fail; docs, Linux
diagnostics and wheel succeed. Download actual source-ee2a7af0 Python 3.13
job 110791892621 once at 10:08:50 UTC: 13,964 passes, five failures, 118 skips,
zero errors and 30 warnings in 482.34s. All 37 inventory and seven record-fixture
cases pass, confirming the preceding correction. The four default raw exports
and the `/proc/self/cgroup` platform probe remain failures. This is not evidence
for the new raw relocation patch, sibling test counts or full CI success.

All relocation documents explicitly retain `redistribution_authorized=false`.
No raw data/native image is committed or publicly uploaded; rights review and
public acquisition/provisioning remain separate requirements. Two exporter
source identities change and a helper is added: prospective release inventories
must include them; historical receipts and archives stay unchanged. Scientific
engine/settings/scores remain frozen. Timing stays deferred without contention
polls/questions, DGX access or unrelated process/service changes. Other-QfO
uncertainty, controlled resources, rights, full release and deposition remain
unfinished; the publication goal is active and incomplete.
