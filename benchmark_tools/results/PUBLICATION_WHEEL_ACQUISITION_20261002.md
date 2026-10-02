# Public Wheel Acquisition and Source Support

Status: 2 October 2026. A copied source component now acquires the ten public
third-party wheels needed by the exact historical inference/reader environments.
It combines them with separately supplied, authenticated OrthoHMM and pip
artifacts and preserves both complete historical locks. This is artifact
preparation, not a complete public-native runtime or new benchmark execution.

## Executed Acquisition

Committed workflow `e15170f366cf669e7784ff3b3600b1bfe12d4fb9` adds
[the wheel acquisition controller](../acquire_publication_wheels.py) and the
explicit `native-wheels` [source profile and commands](../PUBLICATION_SOURCE_COMPONENT.md#historical-inference-and-reader-wheel-acquisition).
The profile extends native preparation with exactly three extra pinned
documents: the admitted 12-wheel inventory and inference/reader locks. Its
seven fixed support files carry metadata, not wheel payloads or raw datasets.
Existing source-only, OrthoBench-input and native-preparation selections remain
unchanged; the current verifier actually accepts the preceding 1,828-file
native-preparation component with its original external index digest.

The [machine execution receipt](publication_wheel_acquisition_20261002.json) is
51,337 bytes, SHA-256
`35c58580b8b74fe7a7636d38ec024c461607ed8f211169ab34d587fb2660b4bd`.
All six actual stages return zero: committed export, copied verification,
public acquisition, ZIP/project-source readback, final copied verification
and historical native-preparation verification.

| Prepared Artifact | Verified Scope |
| --- | --- |
| Public wheels | Ten fresh downloads, 93,684,499 bytes; 20 provider requests including release-specific metadata |
| Supplied wheels | Exact unpublished OrthoHMM 0.5.0 and already acquired pip 26.2.1; neither queried on PyPI here |
| Union | Twelve wheels, 95,645,575 bytes |
| Inference directory | Eleven wheels and unchanged 1,286-byte hash lock, including historical Leiden 0.11.0 |
| Reader directory | Five wheels and unchanged 540-byte hash lock, including admitted Biopython 1.87 |
| ZIP readback | All twelve embedded METADATA name/version identities match the frozen inventory |
| Project source readback | All 33 scientific wheel source members match the separately exported frozen Git blobs |

Wheel names, platform filenames, sizes and SHA-256 come from the fixed admitted
inventory, never the original absolute paths it records. Release-specific
provider metadata must identify exactly that artifact, not a latest version,
different ABI/platform, alternative hash or yanked replacement. Downloads use
the existing bounded allowlisted HTTPS/redirect helper; generated copies/locks
and supplied inputs are rechecked before success. Failures retain partial
files/records without retry, admission or overwriting a previous destination.

Actual execution uses standard-library isolated Python from a fresh copied
source tree outside the checkout. The copied source has 1,833 payloads: 43
scientific and 1,790 workflow files; 1,812 Python files compile without
scientific imports. Its externally retained 741,499-byte index has SHA-256
`c57b58dbb657dc9a6b306a7b6efc1f27a97d282dbdf181373d2a6f1e51345f4b`.
Acquisition's own 27,115-byte `complete.json` has SHA-256
`32186c3ac48b2b9e3a1ea6e9a166924e0b2e5afd3250a7b5256b07645b6fc6b1`.
The complete copied source remains unchanged afterward. The receipt retains
all artifact/metadata/lock/log/checker locations; large wheels are not committed
or publicly redistributed.

## Tests and CI Recovery

The acquisition/base/source/handoff/installer panel passes **260 tests in
21.97s**, zero failures/errors/skips. Its 37,570-byte XML is pinned in the
execution receipt at SHA-256
`b97b27a507438397319a609b6ba1f59b3191bc144a31634a0c00aa62841343c2`.
An earlier 215-case panel overlaps. Fixtures are synthetic; the actual provider
execution and source readback are separate evidence.

While recording this work, preceding source `709b1de6` CI run 37034614556
finishes **failed**. All five test-job logs were downloaded once and read,
with byte/hash retention. Each reports one identical unit failure, 14,460
passes, 119 skips and four declared raw-source deselections; integration was
not reached. Docs, CPU-wheel and Linux-diagnostic jobs report success. Do not
reinterpret the failed matrix or aggregate duplicate test counts.

The failure is a test-only checkout-location assumption: a retained absolute
controller path is split by the current repository root, which is absent in
the macOS checkout. Correction `12f8fdc1` uses `pathlib` relative to the retained
source-root evidence and checks exact ROOT/PYTHON shell assignments. A second
case actually copies the unchanged documents/script into a relocated checkout.
No launcher, runtime pin, environment gate or scientific setting changes.
The [diagnosis and local fix receipt](ci_retained_bootstrap_relocation_20261002.json)
pins all five raw logs and the expanded **371-test** panel, all passing in
22.84s with zero failures/errors/skips. It does not claim the corrected remote
matrix has completed; subsequent source CI must be observed separately.

## Remaining Boundaries

OrthoHMM 0.5.0 is still an unpublished supplied wheel, not publicly acquired
or rebuilt by this controller. Embedded source equality does not prove the
compiled binary's provenance or compatibility on arbitrary hosts. Exact opaque
lock preservation and artifact role sets are not new dependency resolution,
wheel-tag admission or an installed-runtime check. No downloaded code was
imported/executed, no environment installed and no native benchmark rerun.

Conda bootstrap/transitive dependencies, aligner/tree-builder and reader-source
delivery, OS closure, selected-source/rights obligations, complete native-study
handoff, remaining QfO uncertainty, controlled Threadripper timing, final
manuscript and public release/archive/DOI remain open. Reuse admitted full
OrthoBench reproduction and prior installation/fixture evidence rather than
repeating them for unchanged bytes. Same-host relocation is not cross-host or
syscall-sandbox validation. No DGX operation, host-contention poll, new quiet-
window question or unrelated job/service change. The publication goal remains
active, not complete.
