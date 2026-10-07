# Model-Distance Replay V2: Safe Retained Directory Headers

This is a separately versioned repair of the failed first portable replay.
The original source/component/attempt, failed receipt and evidence archive
remain unchanged. Actual v1 rejected the archive before extraction: it contains
183 regular payloads plus 19 zero-sized ancestor directory headers. No selected
feature, score or inference was recomputed in that failed invocation.

V2 permits only zero-sized directory headers that are ancestors of the exact
183 expected files. It rejects duplicate names, unsafe paths, unknown
directories, symlinks, hardlinks and any other member type. Every regular
payload is still byte-checked before extraction; directory metadata is not
restored. All numerical functions and frozen scientific pins remain identical
to v1. The version receipt records the complete four-edit source transformation.

Commit/push and test this source before a new explicitly versioned selected
component and replay. Never reuse the failed v1 namespace or retry it blindly.

```sh
python -B -m benchmark_tools.reproduce_swiss_model_divergence_v2 build \
  --repo . --revision COMMITTED_V2_REVISION --output /fresh/v2-component
python -I -B /relocated/v2-component/replay.py replay /relocated/v2-component \
  --index-sha256 EXTERNALLY_RETAINED_V2_INDEX_SHA256 --output /fresh/v2-replay
```

The component retains the same 12 pinned inputs/readers and three support
files, with the new runner/guide and explicit v2 schema. Before copied code
execution, verify its regular payloads and index against external anchors.
Use Python 3.10 or newer with Biopython 1.87 and NumPy 2.2.6 as declared in
`requirements.txt`; installed 3.10.13 satisfies these requirements locally.

The v1 guide in the source repository records the original attempt's arithmetic
and limits. V2 still invokes only frozen edge/rational functions on restored
inputs, not historical Git/scheduler/binary-dependent admission routines.
All 18 families/563 members/10,765 distances and 54 count records, nine scores,
six contrasts, both TSVs and 15 human table rows must match. No new inference,
scoring, bootstrap draws, uncertainty, defaults, timing repair, independent
confirmation, runtime archive, public clearance, DOI or publication readiness
follows. Actual selected execution remains unproved until a retained receipt
shows it completed.
