# Actual Relocated Model-Distance Arithmetic Replay

The separately versioned v2 source was tested and pushed as
`2ea4be041226bdfe0f06bce5a093b6201f18c16b` before one selected component build
and one actual v2 replay. The original v1 source/component/failed attempt and
the scientific evidence archive remain unchanged.

## Actual Result

- [Build](swiss_model_divergence_portable_build_20261007_v2.json):
  15 flat files, 925,023 payload bytes and an externally retained 4,815-byte
  `REPLAY_INDEX.json`, SHA-256
  `83e3b34b7b99de6f23b6f2c4a4c2cbf28d001c0ac8708b06f427a909346bf3e2`.
- Copied to `/tmp/orthohmm_swiss_model_component_20261007_v2`, with all 15
  payloads byte-checked before execution. Copied runner SHA-256
  `a87e81fe6015e3c12ac29ba2572f1494d3fb2e17e90215d170847b51abc9dbc2`.
- [Actual execution and exact guard command](swiss_model_divergence_portable_execution_20261007_v2.json):
  Python 3.10.13, `-I -B`, Biopython 1.87 and NumPy 2.2.6. Two original-path
  canaries were blocked. No forbidden original scientific source/logical-input
  access or subprocess event occurred during replay. Approved runtime
  site-packages were allowed; this Python-event guard is not OS containment.
- [Complete replay result](swiss_model_divergence_portable_replayed_20261007_v2.json):
  all 183 regular evidence payloads restored and checked; all 18 family
  descriptors, 563 members, 10,765 distances and original median/tie bins matched.
  All 54 count records, nine scores, six contrasts, both TSVs and 15 human-table
  rows checked with the original frozen edge/rational functions.
- Median family cutoff reproduced as `2.24582463005`, with nine families in
  each bin. Approximately 1e-15 summation differences remain within the original
  1e-10 descriptor tolerance and do not change bins or table values.
- [Execution measurement](swiss_model_divergence_portable_replay_20261007_v2.time.txt):
  exit 0, 0.47 seconds, 50,688 KiB peak RSS and zero swaps. These are replay
  costs, not native OrthoHMM timings or isolated-efficiency observations.

## Archived Component

[New portable component](swiss_model_divergence_portable_20261007_v2.tar.gz):
693,861 bytes, SHA-256
`6d995e04580942920301abcbcf3837188c8a8e8480dcd2eb9c13cf30c440fd0c`.
Its 16 regular members are the exact 15 component files plus the externally
anchored index; every serialized payload hash was checked against that index.
Archive modes are normalized to 0644. No directory or link members are present.
The selected replay used the externally checked directory copy, not a second
post-serialization scientific replay. Do not mislabel this as a repeated run.

The component's `README.md`, requirements and copied runner give executable
instructions. Verify the archive, index and runner digests against the external
anchors before copied-code execution; use a fresh destination/output.

## Retained Failure

[Actual v1 failure](swiss_model_divergence_portable_execution_20261007_v1.json)
and [failed receipt](swiss_model_divergence_portable_failed_20261007_v1.json)
remain failed. Its checker rejected the original 202-member evidence archive
before extraction: 183 regular files plus 19 safe zero-sized ancestor directory
headers. The [four-edit source revision](swiss_model_divergence_portable_source_revision_20261007_v2.json)
permits only known ancestor headers, preserving all numerical code/pins and
duplicate/path/type/hash gates. No failed namespace was reused or retried.
The v1 measurement remains exit 1, 0.17 seconds, 30,720 KiB RSS and zero swaps.
Twelve repair tests passed before v2 execution; v1 had 23 prior fixture tests.

## Limits

This is same-function relocated arithmetic replay, not a second independent
implementation, raw-count admission, alignment reproduction, IQ-TREE inference,
model-adequacy validation, ancestral truth or native OrthoHMM reproduction.
Runtime wheels/OS are not archived; another host has not been exercised.
No new benchmark score, bootstrap draw, uncertainty, default, timing repair,
independent confirmation, transitive study closure, public clearance, DOI or
publication readiness follows. The current manuscript/direct-review package
and historical rc5 are unchanged and do not silently acquire this later component.
