# Three Kingdoms Method Input Evidence

## Finding

The retained SonicParanoid2.0.9 run has a different Danio FASTA from the
canonical staged panel. Its native `snapshot.tsv` records SHA-256
`3f43bb141494a64a31bfe8d236fbf6ae54f7112aec4cab83e170a6344c897f7c`,
matching its retained input copy and the raw download. The canonical staged
Danio SHA-256 is
`c0a902d57898468f0df1a77a0773f935ef80d53abc4c3831381736f8d13d4e73`.
The other eleven native snapshot entries match canonical staged hashes.

The [source audit](THREE_KINGDOMS_SOURCE_AUDIT_20260918.md) establishes that
these files differ in seven sequences through removal of29stop markers.
No affected gene is in the scored BUSCO reference, but indirect graph effects
and SonicParanoid's effective internal handling are not established. Do not
claim that the mismatch is harmless or that it explains a score difference.

## Retained Evidence

| Method | Saved input SHA-256 manifest | Retained FASTA copies checked | Result / limitation |
| --- | --- | --- | --- |
| OrthoHMM high sensitivity | None inspected | None in this audit | Metrics record canonical input directory, not per-file historical bytes. |
| OrthoHMM satellite_v2 | All12 staged hashes match | None in this audit | Recorded input identity matches; transformed intermediates not traced here. |
| OrthoFinder3.1.5 full | All12 staged hashes match | 12 | All copied bytes match staged inputs. |
| OrthoFinder3.1.5 sequence-only | All12 staged hashes match | Shared full run | Checkpoint-derived diagnostic, not an independent native inference. |
| ProteinOrtho6.3.6 | None inspected | 12 | All current retained copies match staged inputs. |
| SonicParanoid2.0.9 | Native snapshot | 12 | Eleven staged matches; Danio matches raw, not staged. Snapshot matches every copy. |
| FastOMA0.3.5 | All12 staged hashes match | 12 | All copied bytes match staged inputs, despite .fa filename suffix. |
| OrthoMCL1.4 | All12 staged hashes match | None in this audit | Recorded input identity matches; merged native FASTA not traced here. |

Five ordinary manifests contain the exact12 canonical path/hash entries.
All48 copied FASTAs were hashed; the only mismatch is SonicParanoid's Danio
copy. The native Sonic snapshot also records matching protein counts. Its
reported residue totals are preserved, not independently verified here.
All inspected evidence files and canonical/raw inputs were rehashed afterward.

Saved manifests and current copies are not immutable execution attestations.
The high-sensitivity directory record cannot prove past bytes; conversely,
absence of a hash in that record is not evidence that its input differed.
These distinctions are retained in `three_kingdoms_method_inputs_20260918.json`.
The existing score table is preserved as historical evidence, not promoted
as a fully established uniform-input cohort.

## Selective Rerun Disposition

A fresh SonicParanoid matched-staged-input run is justified to remove the
known mismatch, not to improve its score. Before execution, freeze the same
2.0.9 default-mode command with32threads, all12 exact staged inputs, fresh
copies/cache and the current verified dependency/interpreter environment.
Retain the historical raw-input run under its existing identity. No missing
or unsuccessful method should be silently dropped from comparison.

Freeze conversion and the existing BUSCO-reference scoring rule before
inspection of new outcomes. Validate native completion, all input snapshots,
species-pair outputs, output identifiers and scoring independently before
updating a publication table. Report both input identities and the old/new
results, including unchanged or worse scores. Unless historical dependency
identity is established, do not interpret their difference as the causal
effect of stop-marker removal alone. Shared-host timings remain descriptive.

No new Sonic run was submitted by this audit. Execution and scoring are still
required. The older high-sensitivity input-record gap and other transformed
input paths need separate investigation; this selective rerun would not by
itself prove complete historical input equivalence for all eight methods.

## Validation

```bash
python -m pytest -q tests/unit/test_audit_three_kingdoms_method_inputs.py tests/unit/test_audit_three_kingdoms_sources.py
python benchmark_tools/audit_three_kingdoms_method_inputs.py --repo . --output /tmp/three-kingdoms-method-inputs.json
```

26 focused tests pass; the actual full evidence audit succeeded. Earlier in
this turn, the complete unit suite at4e2e93c passed:4,621tests,9skipped,
86.62seconds. The11 newly added method-input tests were then run in the
focused suite; the older full-suite count does not include them.
