# Current Private Runtime Lookup

Fix a concrete timing-preparation gap, not a timing result. The retained private
OS/helper manifest predates the executor and environment-worker changes.
Both source hashes differ; the lookup inspector itself is unchanged. Replacing
the hardcoded lookup hash in an inventoried executor would invalidate the new
manifest. Preserve the historical pins and add optional, externally pinned
`runtime_lookup` to an explicitly private execution request instead.

The request, readiness review and runtime checker must agree on the lookup
digest. Shared/implicit deployments cannot use the option. Reject malformed,
relative/symlinked or changed records, changed scientific baseline, controller,
command plan, baseline-path coverage, retired roots or private-runtime manifests.
The existing script, 27-run plan, policy-v2 requirement, history/environment/
runtime checks, release gates and no-retry contract are unchanged.

## Executed Checks

At source `b9e07048f6d3b566be2dc3a729b5988986efb7f6`:

- Initial three-module panel: 157 passes in 2.17s.
- After strengthening scientific-binding checks: 164 passes in 2.19s.
- Pre-startup eight-module panel: 319 passes in 8.94s, zero failures/errors/skips.
- Actual private-controller binding rebuild succeeds using the existing tool.
  Its scientific baseline, command plan, controller, baseline paths, retired
  roots and private-tree manifest all match the historical binding.
- Native OS/helper manifest gains 76 Python helper files and changes 21;
  no removed records or non-helper changes occur. Do not overwrite the original.
- Full inventories match before and after the fresh import probe: 39,207
  OS/helper and 18,751 private-runtime records (57,958 total at each boundary).
- One fresh declared-import probe per native interpreter matches its retained
  lookup exactly: OrthoHMM 913 modules, OrthoFinder 1,563. The historical reports
  are reused as the comparison baseline; there are not two fresh probes per tool.

Add eight explicit-lookup executor composition cases after startup: success,
measurement/worker/stream failure, cleanup, selected checker path/digest and
no retry. The same eight-module panel now passes 327 in 9.24s without failures,
errors or skips. Only test bytes change; the executed executor source is still
byte-identical to b9e07048. Preserve the preceding 319-case receipt and source
snapshot separately; no binding rebuild or native probe is repeated.

The controller runs the inventory and inspector; scientific libraries and
their native extensions are imported, but no inference or benchmark scoring
executes. Loader overrides are removed only in copied child environments.
No shared environment, service or unrelated process is changed.

Export the [new lookup receipt](threadripper_private_lookup_explicit_20261002.json)
byte-identically from the actual completed local result:
6,176 bytes, SHA256
`3c244a46786cb0efdf58b15db99b30eb0f166b554e0511d8f07afd20fcc79417`.
The [binding/check/test readback](threadripper_explicit_lookup_validation_20261002.json)
records 20 current evidence identities, manifest deltas and the separately
Git-verified historical test-source identity:
22,634 bytes, SHA256
`3f7334b0f2d584598c3cb665a43278577ffd08fdfe5db6ae76e07bd9e427f37a`.
Latest JUnit: 50,569 bytes, SHA256
`c359da53a2015c85c9d305b027f0735bdec261dbe0e951bc03a7d2a53101142d`.
Retain the exploratory TypeError from treating a metadata records count as a
list; it changed no artifacts. All owned process handles are terminal.

## Remaining Boundaries

No real execution request, passing readiness review, allocation, session,
environmental handoff, observer slowdown or quiet-host proof is created.
Native declared imports are not complete controller import closure, continuous
containment, all workload branches, OS compatibility or security certification.
The standard RuntimeChecker still checks manifests/import lookup around native
work; the request option does not bypass it.

Actual private startup evidence is now current at this source snapshot;
runtime/source drift after this snapshot must still be checked at launch.
The default historical lookup remains unchanged, not newly approved.
27 production and 54 engineering identities remain unstarted. No calibration,
scientific benchmark, main-PDF render or dated archive is repeated. The original
publication goal remains incomplete.
