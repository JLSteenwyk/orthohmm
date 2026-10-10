# Controlled Fragment Native-Stage Relocation

The [actual corrected execution](controlled_fragment_relocation_execution_20261010_v2.json)
completed preparation, relocation, copied-runner replay and independent result
readback. Source `96675b89` was committed and pushed before production.

| Observation | Actual Result |
| --- | --- |
| Pinned component files | 1,783 |
| Pinned payload bytes | 39,428,185 |
| Additional configured original FASTAs | 56 |
| Copied checkpoint FASTAs | 112 |
| Ownership reads / selected native contexts | 352 / 44 |
| Selected cases / reconstructed observations | 53 / 106 |
| Original native checked-input references | 415 |
| Original readback object | Exactly matched |
| Original artifact paths guarded by Python `open` audit | 1,783 |
| Original artifact read attempts during copied replay | 0 |
| Reader runtime | Python 3.12.3; NumPy 2.2.6 |

The [retained manifest](controlled_fragment_relocated_manifest_20261010_v2.json)
is 1,075,354 bytes, SHA-256
`cf42e0e5a9a0486879b976b5cc7848a94cb00d657caacbb64b77a1ba7d04f652`.
The [successful result](controlled_fragment_relocated_readback_20261010_v2.json)
is 14,711 bytes, SHA-256
`0128cd57399dbb4d1bdb7b72660d3be046b098f802059b81eb30162e6d25257f`.
The prepared raw component remains local under
`benchmarks/work/controlled_fragment_relocated_stages_20261010_v2`; its tested
relocated copy and commands are identified in the execution record.

## Failure And Recovery

The [first execution](controlled_fragment_relocation_execution_20261010_v1.json)
failed at ownership lookup. Its plan copied selection/stage pins and checkpoint
FASTAs but omitted 56 configured baseline-original FASTAs. Its source, prepared
component and relocated copy remain unchanged. The prospective `v2` correction
adds ownership pins from retained execution inventories for all 44 selected
contexts. No current file hash was substituted for a missing historical binding.
The original frozen copy writer, independent native reader, NA-table wrapper
and MCL syntax parser are reused, with unchanged kernel bytes. The 74-test
affected/adjacent run passed; the five correction tests also passed after the
final source-freeze assertion was added. Synthetic standalone tests remove the
original input tree before replay. Tests are distinct from this actual run.

## Scope

Historical path labels stay in scientific records. The map binds them to
verified copied files with no original-path fallback; numeric hits, group
partitions, graphs, candidate/event/root reconstruction and exact NA table
serialization remain strictly checked by the unchanged independent reader.
Native kernel bytes and the complete expected result were separately compared.

This closes the specific absolute-path raw-stage diagnostic gap in original
requirement 7.4. It is not new inference, scoring, bootstrap evidence, biological
confirmation, OS isolation, all-tool dependency certification, a replacement
release candidate or proof of full publication readiness. No timing panel was
repeated; shared-host limitations remain. Raw payloads are not committed or
publicly deposited. Follow the [reproduction instructions](../CONTROLLED_FRAGMENT_REPRODUCTION.md)
for the pinned local workflow; summary-only clones do not contain those payloads.
