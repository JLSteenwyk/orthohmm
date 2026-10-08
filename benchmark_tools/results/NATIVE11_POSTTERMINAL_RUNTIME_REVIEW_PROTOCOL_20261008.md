# Native11 Postterminal Runtime Review Protocol

This separately versioned, read-only component addresses the two exact
postterminal OS additions that caused review24033 to correctly reject the old
runtime inventory. It does not rerun native inference or the original full
review, change frozen sources/manifests, remove files, or admit scientific results.

## Evidence Contract

`review_native11_postterminal_runtime.py` binds original request SHA256
`7bf63b80bd5932b9edbd1b2c5ff3fb77f5557e4f6c64e077045d6a50c8d366a1`
and retained classification SHA256
`8718376a789f2175305d83496984ad057aaa559b57469638e97f05ec1fd9537b`.
All classification references, frozen execution sources and original input/session
bindings are rechecked. Fresh accounting must still establish native23985 success.
The original reviewer failures23986 and24033 remain failures.

The unchanged original runtime-review kernel replays before/after inventories,
lookup verdicts, original inputs and copied input bytes. Its explicit checker
hook supplies only the retained, bound first-terminal descriptors. The successor
compares its replayed historical phases to the original partial runtime report;
it does not label the hook output as fresh current identity.

The unchanged inventory kernel then freshly reads all original OS/helper and
private runtime roots. Every original record and all inventory metadata must
remain exact. The OS tree may differ only by the previously pinned regular-file
mode/size/SHA256 descriptors for `/usr/bin/lftp` and `/usr/bin/lftpget`; the
private runtime must match entirely. The retained accounting and package-log
line are parsed to reproduce the strictly postterminal chronology. The original
reviewer's end is an upper bound, not an invented precise inventory-check time.

The new schema `native11_postterminal_runtime_component_v1` explicitly reports
current equality to the old OS snapshot as false. It keeps continuous-integrity,
ordinary-full-review success, scientific admission, accuracy, resources,
next-identity authorization, retry and publication readiness false. A separate
composed full review must still satisfy resource/environment/output contracts;
ordinary frozen conversion/history consumers cannot accept this component.

## Execution

Commit and push source, tests and this protocol before one execution. Use the
existing scientific Python3.10 and sanitized environment, one native thread per
library. Invoke `-B -m benchmark_tools.review_native11_postterminal_runtime
--source-sha256 <committed-source-digest>`. The fixed fresh destination is
`benchmarks/work/native11_postterminal_runtime_review_20261008_v1`.
No new Slurm inference allocation or native repeat is required. Retain ordinary
component failures and diagnose them; do not automatically retry. The program
refuses an existing destination and prints only its resulting record.

## Tests

Focused fixtures exercise the actual unchanged bracket/lookup kernel and real
inventory walker. They reject changed/missing originals, unexpected/altered
additions, private differences, metadata drift, duplicate inventory paths,
altered phases/descriptors, wrong lookup and input bytes, invalid chronology,
unsupported package evidence and duplicate accounting. They check historical
versus current labeling and ordinary conversion rejection. These tests do not
claim that the production component or composed full review has executed.
