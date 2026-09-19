# Held-Session Root Context Panel 22020

## Result

The separately identified full follow-up completed `0:0` in 4 minutes
26 seconds. All 12 workloads and their raw observations pass the v2 audit,
including service command, CPU dose, overlap, scope, exit/removal, fixed
order, source/runtime checks and scheduler allocation. Fresh-directory
archive extraction and replay also validate all 12 trials. Panel 22019
remains separately retained as incomplete; no missing repetition was replaced.

| Condition | Valid trials | All/common intervals | Original flags, all/common | Narrow flags, all/common |
| --- | --- | --- | --- | --- |
| Idle | 3 | 63 / 57 | 0 / 0 | 0 / 0 |
| Steady | 3 | 63 / 57 | 0 / 0 | 0 / 0 |
| Churn | 3 | 63 / 57 | 0 / 0 | 0 / 0 |
| User-contended | 3 | 63 / 57 | 60 / 57 | 60 / 57 |

User-slice CPU over the three contended enclosing observation spans was
18.918252, 19.119470 and 19.038010 seconds. All exceed the unchanged
five-CPU-second response criterion. Across the nine idle/native-only trials,
the same aggregate ranged from 0.013090 to 0.021242 seconds. The intended
response was therefore observed in each valid positive control, but this
does not identify every CPU source or recover isolated task CPU exactly.

The description retains every interval, its actual enclosing read durations,
signed residuals, separate host categories, root membership changes and
original/narrow flags. Fixed-block median differences are descriptive only;
dependent intervals are not independent replicates for inference.

## Session And Environment

Source commit `4c57182`; recipe committed in `bbddb09`. All 512 deployed
files matched source archive bytes and sizes, with 515 recipe records.
Recipe SHA-256:
`a0a2fdbb98fe75b057b2acfb3d9eb94ac1bb6acab903c5a5f71a30ce32a4b627`.
Fresh deployment/output names are `root_context_controls_recipe_v2` and
`root_context_controls_v2` under the established DGX project root.

The bounded launcher retained an empty scheduler queue and used one waiting
SSH session as specified in `ROOT_CONTEXT_SESSION_PROTOCOL_20260919.md`.
Its receipt independently matches job 22020 and brackets the terminal Slurm
lifetime, using the recorded America/New_York scheduler timezone and allowing
for one-second printed timestamp resolution. The scheduler recorder completed
46 polls with zero observation errors. No extra SSH inspection or transfers
were made during the panel. No lingering setting was enabled.

The bounded historical journal query shows one manager start at the beginning
of the job and no manager stop/restart within the retrieved interval. It
also records 50 scheduled restarts of the unrelated
`samwise-daemon-samwise.service`, which fails because its working directory
is missing. This activity was not stopped or excluded. Together with the
held session/client, it is part of the observed environment, not evidence
of a background-free machine. User-slice responses are aggregate values.

This validates the bounded monitor response under these conditions. It does
not identify the cause of prior native-tool flags, establish general process
attribution, quantify added context-collector overhead on native tools,
prove non-CPU isolation or admit scientific comparative timings.

## Evidence And Reproduction

- Raw/source archive: `root_context_archive_22020.tar.gz`, 1,344 files,
  SHA-256 `936d437fcb2d08a35d17ba82a8a48ed18352ce0b00f16dea6bb8be50f523be3f`.
- Submission and scheduler receipts: `root_context_session_receipts_22020.tar.gz`,
  SHA-256 `1499f5443e36e5ba025c9cdf0811b6defda6742e9361139bd2ea24c1b5164e2f`.
- Full audit: `root_context_audit_22020_20260919.json.gz`, SHA-256
  `9b850ca512ac8988f8cdd8b246bace521cc611d842b1e1da7a52598e2726c1db`.
- Description: `root_context_description_22020_20260919.json`, SHA-256
  `4edad9736fe154c7d6b3face89b8fecadcc44578ab86b1afc7b0a03fa9e0e5af`.
- Receipt audit: `root_context_session_receipt_audit_22020.json`, SHA-256
  `d3304832184b3ed77ad3086010a01f006fede98bbfbba030eab629774cc56358`.
- Journal: `root_context_user_manager_journal_22020.json`, SHA-256
  `d1bd0f3495491a17d63954c370f1120dc48d4fb02bb0851dd1d9617377ebdd42`.

Extract the archives into fresh directories. Run
`benchmark_tools.audit_root_context_controls` with `--deployment v2`,
`--job 22020`, the retained recipe/hash, terminal scheduler record and extracted
raw archive. Separately run `audit_root_context_session_submission` against
the extracted receipts, recipe and scheduler, with
`--scheduler-timezone America/New_York`. Run `describe_root_context_controls`
with the full audit hash above. Generated evidence paths record each replay
location and need not equal the original machine's absolute paths.

Validation: 488 focused tests passed before deployment. The broader unit
suite passed 7,688 tests with nine skips in 184.93 seconds while the panel
ran on the separate host. After adding real-result regression fixtures,
49 targeted receipt/audit/description tests passed. No accuracy defaults or
queued corrected-QfO workflows were changed.
