# Corrected FastOMA Configuration And Resource Probe

## Outcome

The installed Nextflow 22.10.8 and pinned FastOMA Docker image successfully ran
a one-task cgroup-v2 probe on bizon. The final wrapper contains one immutable
image reference, `--cpus 1.0`, `--memory 256m` and `--network none`. Inside the
container, `cpu.max` was `100000 100000` and `memory.max` was `268435456`.
The trace records exactly one completed, uncached task with exit zero.

`fastoma_corrected_resource_probe_20260918.json` records native artifacts,
effective configuration, validated selected settings and image identity.
SHA-256: `f2f2e8d62acfee09311536ceb458746de8824936f4a10f07640b1753c7c5e761`.
Configuration SHA-256:
`e16e87d5687bc6e88fc793626418effc34a34271261b67a36721c32cfc605a95`.
Forty-one focused tests passed, including resource-limit mutations, mutable
images, duplicate flags, changed algorithm settings and incomplete traces.

This is not biological inference, a corrected accuracy result, a full runtime
freeze or proof of aggregate resource enforcement. Artifact hashes were
recorded after the observed probe. Production admission remains separate.

## Configuration

`fastoma_corrected_execution.config` explicitly includes the retained host
workflow configuration and historical QfO override. It retains 180 maximum
task CPUs and a 700 GiB task-memory pool, OMAmer 24 GiB/28 concurrent tasks,
UniProt IDs, forced native ortholog-pair generation and the local pinned LUCA
database. The LUCA symlink resolves to the retained database, whose complete
9,375,445,304-byte content was rehashed against the asset manifest.

Filtering remains `col-row-threshold`, row/column gap ratios 0.3/0.5 and five
representatives per HOG. The configured minimum sequence length remains 40;
this does not add a missing command-line option to the retained host workflow's
`infer_roothogs` process. That script remains the historical version, not the
different copy inside the image.

Collection and pair extraction now request 16 CPUs/280 GiB each, matching the
successful historical recovery allocation rather than repeating its failed
initial low-memory collection. This is a prospective resource adjustment,
not an algorithm change or a claim of runtime equivalence. Production should
reserve 720 GiB total, including 20 GiB controller/overhead allowance above the
700 GiB task pool. The allowance is not a measured upper bound.

The configuration pins the immutable Docker digest and enables a task trace.
Use Nextflow `-C` without `-profile docker` to avoid merging unrecorded user or
launch-directory configuration or selecting the mutable image-tag profile.
Nextflow documents the fixed-configuration behavior of
[`-C`](https://www.nextflow.io/docs/latest/config.html). Its current
[local-executor documentation](https://docs.seqera.io/nextflow/executor/local)
distinguishes scheduling requests from runtime enforcement; the retained probe
tests the actual installed 22.10.8/Docker behavior rather than assuming current
documentation applies identically to that release.

## Retained Attempts

- `benchmarks/work/fastoma_corrected_resource_probe_20260918/`: failed before
  task submission. The tiny standalone script inherited a time-limit closure
  requiring `check_max`, a function provided by FastOMA's full workflow.
- `fastoma_corrected_resource_probe_v2_20260918/`: succeeded after setting an
  explicit one-minute probe limit. Its wrapper revealed that Nextflow already
  emitted a hard CPU limit, making the added CPU option redundant.
- `fastoma_corrected_resource_probe_v3_20260918/`: final successful fresh run
  after removing the redundant option. No prior work directory was resumed.

Warnings about unmatched biological process selectors are expected for this
single-process probe. It does not instantiate or validate all dynamic memory
closures from the complete biological workflow. Raw failure evidence remains
retained; no unrelated container/job was stopped and no DGX scan was used.

## Remaining Gates

Require the admitted corrected OrthoFinder tree, fresh proteome/tree copies,
full Java/Docker/Nextflow runtime identity, an exact launch manifest and
scheduled resource verification before inference. The configuration still has
machine-specific absolute paths; portable release packaging is unfinished.
Native completion, pair integrity, conversion and six-endpoint scoring need
their own subsequent checks. No corrected FastOMA inference job is queued.
