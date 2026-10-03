# Portable Resource Reporting Component

The [versioned local archive](shared_resource_reporting_component_20261003_v1.tar.gz)
contains 15 payload files and one externally pinned manifest. It preserves the
three-attempt reviewed snapshot, attempt/cell tables, original figure outputs,
project license, reporting requirements and unchanged reporting source.
The [result receipt](shared_resource_reporting_component_result_20261003.json)
records actual fresh extraction, standalone verification and copied-reader
replay. This is a reporting component, not native measurement reproduction or
a complete 27-attempt study release.

Archive: 315,670 bytes, SHA-256
`d14fc99f9ee4580e4b3a57395420f52148b825b2cb2cf0ba7716b57a11e8f6bc`.
Manifest: 2,716 bytes, SHA-256
`e3d73bcb05661b84306114e01cf8aa51b15d41c58025ecf9313c46cc4f0c39c4`.
Sources are committed at `4ac940d4`; the reader selects existing pure functions
and constants using the repository's established AST-selection pattern. It
does not execute the original collectors, source imports or scheduling APIs.

After checking the archive hash, extract into a fresh directory and use its
own reader. Standard-library verification needs no project checkout or
scientific packages. Plot replay uses the recorded Matplotlib 3.10.8 reporting
environment; this run also used NumPy 2.2.6. An existing compatible environment
is sufficient; the archive is not a full runtime or dependency wheelhouse.

```sh
python -B -S component.py verify --directory /absolute/restored/component \
  --manifest-sha256 e3d73bcb05661b84306114e01cf8aa51b15d41c58025ecf9313c46cc4f0c39c4
python -B component.py replay --directory /absolute/restored/component \
  --manifest-sha256 e3d73bcb05661b84306114e01cf8aa51b15d41c58025ecf9313c46cc4f0c39c4 \
  --output /absolute/fresh/reporting-replay
```

All three table files reproduce byte-for-byte and the PNG pixels match exactly.
The two eligible observations, excluded first attempt and missing repeats are
unchanged; no incomplete median/range is fabricated. Replayed PDF/SVG outputs
are generated, but their binary identities are not asserted equal because of
format metadata. The original outputs are also preserved in the archive.

The copied replay rejected an original-checkout canary, recorded no later
original-project open events or external process execution and imported no
project modules. The reporting interpreter's installed prefix was explicitly
allowed. These Python audit events do not establish OS containment, complete
native-I/O observation, cross-host portability or an independently provisioned
runtime. Original absolute evidence references remain provenance only and are
not opened to revalidate cgroup counters, runtime or native outputs.

The initial 16 focused tests pass, including actual table/pixel replay and
corruption, inventory, scope and overwrite refusal. The retained-archive
regression additionally binds the actual payload hashes and scope. The combined
component, table and figure suite passes **73 tests in 4.07s** on resumption.
No inference,
score, uncertainty result, biological validation, isolation or speedup claim
is introduced. Full resource/figure/manuscript/archive reconciliation remains
pending; this Git-managed component is not a public archival deposition or DOI.
