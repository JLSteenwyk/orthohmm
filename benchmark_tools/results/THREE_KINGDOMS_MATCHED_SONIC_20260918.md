# Matched-Input Three Kingdoms SonicParanoid Rerun

## Frozen Design

The historical SonicParanoid run used raw Danio sequences, whereas the
canonical staged panel removes 29 stop markers from seven sequences. The
historical result remains retained and labeled; this is a fresh matched-input
run, not a pure causal test of that transformation. Current dependency
identities do not establish the historical April runtime.

The frozen plan is `three_kingdoms_sonic_matched_commands_20260918.json`,
SHA-256 `5c3ce4608b54e8e550f0ea855754879211b842f0175b2e8e864cf4b5c2998283`.
It binds all 12 staged FASTAs (443,217 proteins), audited source records,
SonicParanoid runtime and tool inventories, and the existing conversion and
BUSCO scoring scripts. No old searches are reused. Native default mode is
run with 32 CPUs, 192 GiB, and a 72-hour limit on bizon.

The output root is `benchmarks/results/three_kingdoms_sonic_matched_v1`.
The launcher refuses an existing root, checks copied input bytes, records
native exit status and GNU time output, and rechecks runtime/input hashes.
Preflight completed successfully without creating the inference root.
The 36 focused launcher/runtime tests passed.

## Admission Requirements

Native success is not accuracy admission. Before updating comparative tables:

1. Verify terminal scheduler state, actual resource allocation, native input
   snapshot and output identities against the frozen plan.
2. Validate native orthogroup identifiers, membership, and input ownership.
3. Use the existing SonicParanoid group-co-membership normalization, not the
   native pair-table semantics used for QfO.
4. Run the unchanged BUSCO scorer against the same 255-group reference and
   independently verify pair-count arithmetic and denominators.
5. Report the new result alongside the historical mismatched-input result,
   with runtime-version and reference-scope caveats.

Shared-host time and memory are descriptive, not dedicated scaling evidence.
Scheduling after QfO numeric admission is a resource-priority decision, not
a scientific dependency. A failed QfO admission does not invalidate this
independently frozen Three Kingdoms input panel.
