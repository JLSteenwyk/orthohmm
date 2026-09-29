# Prospective Private Deployment Amendment

Prepared a [new baseline](threadripper_private_baseline_20260928.json) and
[new command plan](threadripper_private_commands_20260928.json), retaining the
originals. The [amendment receipt](threadripper_private_deployment_20260928.json)
links both originals, the patched installation and successful native fixture
parity. This is deployment preparation, not runtime binding or execution approval.

Only OrthoHMM's interpreter entrypoint record and package environment metadata
change in the baseline. Only `native_argv[0]` and `configuration.argv[0]` change
in each of the 18 OrthoHMM run rows. The nine full OrthoFinder commands are
unchanged. All 27 run identities, datasets, source hashes, settings, output paths,
resource limits and rotated method order remain unchanged. Historical
`original_native_argv` fields intentionally retain their provenance values.
The patched installation's live isolated inventory still matches all 30 pins.

## Existing Verification

The unchanged `verify_environment` implementation passes from the frozen
native working directory with the baseline environment overrides and a fresh
verification-only Numba cache. The [receipt](threadripper_private_verification_20260928.json)
records that context. This exercises existing native binary/profile checks,
frozen source/adapter/tool records, core revision, and both interpreters'
package inventories. It does not establish transitive runtime-tree coverage.

An [earlier context error](threadripper_private_verification_wrong_cwd_20260928.json)
is retained: calling that verifier from the mutable checkout exposed local
OrthoHMM 0.5.0 distribution metadata through the current-directory search path,
causing the expected package-inventory rejection. No package changed. Using
the plan's frozen cwd removes that extra metadata. Future controller launches
must set the documented cwd and environment before verification.

An initial amendment-generation check also assumed fixed method blocks; it
rejected the actual rotated plan before writing files or launching work.
The corrected check explicitly validates the frozen rotation by repeat and
dataset size. Six unit tests pass, including identity/interpreter rejection
and proof that reversing the two intended argument edits reproduces the entire
original plan without mutating it.

## Next Integration

Inventory the private interpreter/base and collector runtime, establish fresh
repeated native lookup evidence, and bind these new baseline/plan bytes to the
collector v5 fixture workflow. Do not reuse the old shared-prefix binding or
claim this amendment has already passed those checks. The controller still
needs an explicit dependency scope; current preparation used shared Python.
Final resource accounting and verified whole-run quiet-host eligibility remain
separate requirements. No production identity, native inference rerun or DGX
access occurred in this amendment turn.
