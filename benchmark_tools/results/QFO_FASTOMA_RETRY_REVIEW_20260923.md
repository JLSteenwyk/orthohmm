# Corrected QfO FastOMA Retry Review

This is a read-only diagnostic review, not native-output admission or an
accuracy result. Original task artifacts and failed admission 21741 remain
unchanged. No failed trace rows were removed and no validator was weakened.

## Attempts and Commands

The trace contains 358 COMPLETED/0 attempts and two FAILED/1 attempts.
Nextflow's log explicitly records the two corresponding resubmissions.

| Process display name | Failed task/hash | Successful retry task/hash | Actual batch | CPU limit | Memory limit |
| --- | --- | --- | --- | --- | --- |
| hog_rest (137) | 258 / 88/b76e7a | 356 / 86/252823 | rhogs_rest/220 | 1 in both | 12 then 24 GiB |
| hog_big (7) | 88 / ac/a2d49a | 357 / fe/c2a4b8 | rhogs_big/14 | 4 in both | 12 then 24 GiB |

Within each pair, `.command.sh` is byte-identical (`diff -u` exit zero).
Both attempts use the same absolute batch path, checked species-tree name,
output name, filtering parameters, and parallelism setting. All four
wrappers pass the existing pinned-image, no-network, resource-limit check.
This is command identity, not proof of historical input-byte identity.

The failed regular task mounts only the `work/88` subtree into Docker;
its input batch is under `work/b6`, outside that mount. Its stderr reports
FileNotFoundError for precisely that batch. The retry mounts the whole work
directory and succeeds. This is strong evidence of an input-visibility
failure, not a biological inference failure. No claim is made about how
Nextflow selected the original common-prefix mount or whether other
unrecorded system circumstances contributed.

The large task reports a killed MAFFT subprocess. Its retry doubles the
memory limit and succeeds, but a killed subprocess plus this resource
change is not sufficient evidence to label the failure OOM. Memory limits
are not measured peak RSS. The native run's elapsed time includes retries;
neither attempt should be removed from failure or resource reporting.

## Coverage and Collection

Read the complete TSV with `csv.DictReader`, select COMPLETED/0 rows, and
resolve task directories with `audit_fastoma_tasks.task_directory`:

- Successful process counts: check_input 1, omamer_run 78, infer_roothogs 1,
  batch_roothogs 1, hog_big 40, hog_rest 234, collect_subhogs 1,
  extract_pairwise_ortholog_relations 1, fastoma_report 1.
- Parse each successful HOG `.command.sh` with `shlex.split(comments=True)`
  and extract `--input-rhog-folder` using the existing `option` helper.
  The resulting unique absolute-path sets match exactly the 40 children
  of the batching task's `rhogs_big` and 234 children of `rhogs_rest`.
- The collector's `pickle_folders` contains 274 entries, all symlinks.
  Resolving them strictly yields 274 unique directories, exactly the
  `pickle_hogs` directories from successful HOG attempts. Neither failed
  attempt's output directory is included.

Batch task directory: `work/b6/9e7f4935d8648584fa16213662b9dc`.
Collector directory: `work/a0/fa6657994a22c450ebd85082309b82`.
Paths here are relative to `benchmarks/results/qfo_corrected_fastoma_v1`.
These checks inspect paths and coverage, not pickle contents; no pickle
was deserialized. They do not establish final output biological validity.

## Evidence Checksums

SHA-256 of retained evidence:

| Artifact | SHA-256 |
| --- | --- |
| run/trace.txt | 13def3f70ccadbf806d7f385c3f34c4de2e218eefa59ece945ed071ce5c52c38 |
| run/.nextflow.log | 8286f67e9465688ed63091bf08e903d80d350ff4afd0d3a0b16a6e129d9cad45 |
| Both regular-task command scripts | b5a88182352c3a0966b339a14f38ff336bb4a6aac2ebf1899036ac44a69bfb5f |
| Both large-task command scripts | 11adbf5d7dd8c989eb0222047c4aedcd6e235ae79da90f133407d09f329a6bde |
| Failed regular wrapper | 50b8f1f35c63cabf3bc48376493b74a2709360e6604b4614583b8b8868a985c4 |
| Successful regular wrapper | b55203864b247db3e1ed39fff0f281c749a97adb90fa9514f7ba37007ba6ac98 |
| Failed large wrapper | 2f19eff62af9643f20bdd65846a1793dbc628e7c2e492c93454d67b73aeeccf7 |
| Successful large wrapper | 061577c6d33fbaceb418abb367b559e4e5dd8451f3bb1e8613d7b98667b7e803 |
| Collector command script | e0e9759a24654a8c3878e43b34b59cf1940c2d6a26be7e9bdcbbb9310807d857 |
| Collector wrapper | 48f536838c2356160172ef7fd660452210ffbec7c63c12e0ff1c4576b048b8fc |

## Remaining Admission Work

Implement and test an explicit retry-aware audit that retains every
attempt, binds retry command/input identities, validates absolute batch
paths against the batching output, and checks exact successful collection.
Keep rejection of arbitrary failed/cached/duplicate tasks. Verify the
frozen retry policy, staged inputs, and published outputs using the
existing independent admission requirements. Do not treat this report as
permission to skip those checks or rerun the expensive native workflow.
