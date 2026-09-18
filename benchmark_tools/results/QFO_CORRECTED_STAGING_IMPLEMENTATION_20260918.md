# Corrected QfO Staging Prepared

`stage_qfo_corrected_inputs.py` implements the separately named canonical
staging step in `QFO_CORRECTED_RELEASE_PROTOCOL_20260917.md`. It has not been
executed against the production archive. Acquisition 21687 remains live;
compatibility audit jobs 21688/21689 are still pending. Passing fixture tests
does not establish empirical release compatibility.

## Preconditions And Behavior

The caller must supply exact SHA-256 values for both completed audit reports
after reviewing scheduler completion, pinned executors, logs and provenance.
The staging command does not perform those scheduler/source-code reviews on
the caller's behalf. It rechecks all report input file identities, original
manifest/mapping/alias/native-database pins, archive size and common source
identities. The two reports must agree on every canonical file and coverage
statistic, including gzip EOF validation.

Acceptance requires exactly 78 expected proteomes, changes confined to the
Xenopus canonical FASTA, 984,137 uniquely mapped sequences with complete
native reference coverage, all 14 formerly missing SwissTrees accessions
recovered, and no unexplained sequence difference. The native classification
is recomputed from its difference records and compared with the reported
classification. B/O/U/Z-to-X representation matches remain separate from
exact matches. The command does not alter residues or silently relax gates.

Extraction creates a new directory exclusively. It streams only allowlisted
regular canonical members to flat filenames, checks exact member paths,
sizes and SHA-256 values, rejects duplicate/unexpected members and unsafe
paths, and drains gzip through EOF. It does not use unrestricted archive
extraction. DNA/additional files are not canonical inputs. All report/source
and staged file identities are rechecked before the success manifest is
written. A failure leaves a partial directory without a success manifest;
the command will not overwrite that directory on another invocation.

The success manifest explicitly says
`corrected_inputs_staged_pending_inference_freeze` and
`inference_authorized: false`. Independent staged inventory and an immutable
executable experiment manifest are still required. No inference, score
replacement or corrected accuracy claim is produced by this tool.

## Invocation After Review

```sh
python benchmark_tools/stage_qfo_corrected_inputs.py \
  --comparison benchmark_tools/results/qfo_corrected_archive_comparison_20260917.json \
  --comparison-sha256 "${REVIEWED_COMPARISON_SHA256:?}" \
  --sequences benchmark_tools/results/qfo_corrected_archive_sequences_20260917.json \
  --sequences-sha256 "${REVIEWED_SEQUENCES_SHA256:?}" \
  --destination "${NEW_CORRECTED_INPUT_DIRECTORY:?}"
```

Obtain the required environment values only after reviewing the completed
empirical reports and choosing a new destination. The shell refuses unset or
empty values. `--help` is verified. No production invocation is queued.

## Tests And Concurrent Jobs

36 tests pass across staging, archive comparison and archive/native sequence
comparison. They cover compatible reports, representation-only vs unexplained
differences, missing coverage, unexpected change scope, conflicting reports,
missing SwissTrees accessions, source-identity disagreement, exact unchanged
extraction, no overwrites, archive symlinks, path traversal, absolute paths,
duplicate/missing members, checksum failure and truncated gzip.

At this milestone, archive acquisition reached 1,444,024,320 bytes out of
2,648,666,198. Original QfO reconciliation 21671_1 is live. Dedicated DGX task
21656_8 completed with scheduler exit 0:0 in 1:36:10; task 21656_9 is running.
That is nine scheduler-completed timing tasks, not nine fully admitted
resource measurements. No unrelated jobs were stopped.
