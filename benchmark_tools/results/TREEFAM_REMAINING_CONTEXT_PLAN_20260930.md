# Remaining Public QfO Context Targets

The [six previously inspected tags](TREEFAM_QFO_DEVELOPER_CONTAINERS_20260930.md)
remain retained and will not be queried again. The complete saved public tag
listing contains 40 names. Inspect the other 34, including 2020.1.1-2020.1.3
and historical commit/developer aliases, under the same bounded copied-context
hypothesis: original inputs absent from Git might have entered a Docker build
context. Tags alone do not establish data-release or build provenance.

The [prospective plan](treefam_remaining_context_search_plan_20260930.json)
fixes the target list and source bytes before inspection. Fetch public
manifest/config records, identify COPY/ADD-to-/benchmark layers and verify
every selected blob's size and SHA256 against its descriptor. Reuse matching
raw blobs from the two prior download directories; corrupt or indirect cache
files fail, rather than silently triggering a replacement download. Deduplicate
layers shared by aliases. Preserve all tag outcomes and unresolved cases.

No container execution, archive extraction, contacts or inferred family labels.
Keep the 32 MB per-context and 96 MB new-context download limits. A registry
429 stops further attempts, without retry or bypass. Oversized/missing context,
unsupported platform, malformed history or transport failure is unresolved,
not evidence that original files are absent. Listing/token failures may stop
the operation before a report; retain acquired files and do not restart merely
because observation ends. Other layers and registries remain outside scope.

The helper now supports explicit tags and descriptor-verified cache reuse;
18 focused tests pass, including actual mocked registry-to-inventory flows,
positive NHX detection, cache reuse/corruption, safe names, timeout retention
and rate-limit stop. This changes only the reference-search helper, not the
frozen OrthoHMM method, scores, native commands or timing identities. Existing
scientific/runtime/calibration receipts remain historical; no current timing
recipe or execution-readiness approval is implied by the new source.

Commit/push the helper, tests and fixed plan before one new inspection. On
completion, independently check raw records, descriptors and complete tar
inventories before drawing bounded conclusions. The original full TreeFam-A
trees/mapping are still missing at this prospective checkpoint. Controlled
timing remains deferred, with no new host probe or coordination question.
