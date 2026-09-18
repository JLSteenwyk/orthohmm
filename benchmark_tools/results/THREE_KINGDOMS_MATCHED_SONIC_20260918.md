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

Submitted as Slurm job **21795** from frozen executor revision
`7c3d0784d1a22da0a0e37db6860977c6649fd2ca` in
`benchmarks/work/publication_three_kingdoms_sonic_matched_v1`.
Controller inspection confirms PENDING, `afterany:21792` unfulfilled,
32 CPUs, 192 GiB, three days, bizon, no requeue and zero restarts.
Controller SubmitTime is `2026-09-18T13:10:52` (timezone not inferred).
No native output or accuracy result exists at submission.

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

## Native Group Validator

`benchmark_tools/validate_sonicparanoid_groups.py` now validates exact full
FASTA identifiers, species ownership, snapshot input hashes/protein counts,
unique native group IDs and membership, declared group/species counts and
seed-count bounds, and bijective membership equivalence with normalized
groups. It rejects unexpected columns and truncated/extra-field rows.
Files are hashed before and after validation. The validator does not by
itself validate scheduler state, runtime identity, or score arithmetic.

The historical retained input copies, snapshot, native table and normalized
groups pass this validation: 12 species, 443217 input proteins, 19853 groups,
288562 grouped proteins, and 154655 proteins outside the selected group
table. Report: `three_kingdoms_sonic_native_groups_20260918.json`.
Full pipe-delimited Xenopus identifiers are preserved without accession
stripping. This uses historical input copies, including raw Danio; it does
not establish equivalence to staged inputs or admit a matched-input score.
Nineteen focused validator/normalizer tests pass. The new job still requires
its own terminal, provenance, conversion and scoring validation.
