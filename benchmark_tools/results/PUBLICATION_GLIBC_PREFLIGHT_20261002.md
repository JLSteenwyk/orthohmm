# Assembled Glibc Preflight

Status: 2 October 2026. Executed preflight-only checks, not installation,
new native integration, scientific evaluation or controlled timing.

## Enforced Constraint

The [declared ABI inventory](PUBLICATION_NATIVE_ABI_20261002.md) identified
`GLIBC_2.34` requirements in all 32 compiled MAFFT helpers. Current assembled
execution now requires `--abi-inventory` and its external SHA-256, tied to the
same externally pinned assembly. Before creating output or installing packages,
it checks the running controller's numeric GNU libc version, Linux/x86-64,
declared interpreter presence and absence of controller loader overrides.
Legacy execution without an assembly flag remains unchanged; mixing the new
ABI arguments into that route is rejected.

The inventory's aggregate requirements are used conservatively, including
weak versions. Unknown/nonnumeric GLIBC names fail rather than being declared
satisfied. Observed controller glibc is obtained through the operating-system
[confstr interface](https://docs.python.org/3/library/os.html#os.confstr).
This is only a numeric floor check: it does not establish which providers the
private base/assets load, CPU compatibility, static/dlopen dependencies,
non-glibc version satisfaction, security or complete OS closure. Loader file
presence and hashing are not loader-resolution evidence.

`--preflight-only` checks the existing input, asset, lock and base-binary
bindings, rechecks them and the raw inputs before finalizing, and writes only
`preflight.json`. It never launches the base package probe, installations,
inference or scoring. It is not a reusable permit: normal assembled execution
performs fresh pre/post checks and retains the other runtime/scientific gates.

## Executed Checks

Initial source commit `89caa096ea2c36ea2027140385be0d8a5d5a7ea0` introduces
the guard/route. The first interactive preflight refuses inherited CUDA
`LD_LIBRARY_PATH` before creating output or launching stages. Preserve this
early refusal as a retrospective terminal-tool transcription, not a fabricated
captured log or independent process-accounting receipt. The host is unchanged.

A separate child-only environment removes `LD_LIBRARY_PATH`, `LD_PRELOAD` and
`LD_AUDIT`. The preflight then succeeds on the same relocated assembly and
retained 16-gene fixture, using the exact private-base binary pin. It reports
controller `glibc 2.39`, declared floor `2.34` and 119 watched file bindings.
Only the preflight JSON is created; no inference/environment directory exists.
Current base-package/runtime probing is explicitly unexecuted.

Review identifies a missing final raw-input recheck in preflight-only mode.
Commit `c2bc6ebdc3d5e2637a7190acc097bd34a348ecb1` adds it and a changed-input
rejection test. One further preflight executes that concrete new check with
the same unchanged data/assembly and clean child environment. The first
successful receipt stays unchanged. There is no native or statistical rerun.

The input-recheck five-module panel passes 128 tests in 2.43s, zero failures/errors/skips.
Earlier overlapping 124/127-case panels pass and are not added to this count.
Tests cover old/unknown GNU libc, unknown version names, platform/assembly/
digest mismatch, loader overrides/unavailability, changed evidence, changed
end-of-run host checks, mandatory route arguments, and no-stage preflight mode.
The final 18,135-byte JUnit receipt is
`benchmarks/work/publication_glibc_guard_final_tests_20261002.xml`, SHA-256
`88f3b1b36f192fd353da2280cf6b6ae4cabab9b617f983950611ec3e72310805`.

Readback verifies 132 current file identities, all seven raw fixture inputs,
both logs/preflight reports, source snapshots and bounded flags. The older
controller record belongs to its recorded Git revision; it is not claimed to
match current mutable bytes after the new input check. No new host observation
or preflight invocation is performed by readback.

A final CLI-only correction at `bf915f2db5ca18216f41da2ec5df684c9d56592a`
rejects `--preflight-only`/ABI arguments combined with the internal scoring
worker, instead of silently ignoring them. The latest panel passes 131 tests
in 2.40s, zero failures/errors/skips, including three actual subprocess
argument-rejection cases. Its 18,597-byte JUnit receipt is
`benchmarks/work/publication_glibc_guard_cli_tests_20261002.xml`, SHA-256
`8151c85a47da7531bc39e2485931925ef3bba80824abb038dc48bedeaebe2508`.
Separate binding readback proves all other controller bytes and the pure guard
unchanged since the recorded preflight, with 118 other watched files and all
seven raw inputs still matching. The 132-identity readback remains dated
before this CLI amendment, not a claim that its controller digest is current.
No host query/preflight is repeated just for argument-conflict rejection.

| Evidence | Bytes | SHA-256 |
| --- | ---: | --- |
| [First clean-environment preflight](publication_glibc_preflight_20261002.json) | 7,722 | `01c502f152d2a444f71b8c3bfb525d74077472dc411fe2004cc946ec15e08903` |
| [Final-input-checked preflight](publication_glibc_preflight_final_20261002.json) | 7,274 | `75bf519f0de3947d3888c546d0117c2a76fe8295af6eba8117c57039524e622d` |
| [Readback](publication_glibc_preflight_validation_20261002.json) | 2,709 | `4b258b672bfb35a0884fe1ef3a725c25a1ec8abb647c7e9e3ac598bbb90ed40b` |
| [CLI-only correction binding](publication_glibc_cli_gate_validation_20261002.json) | 3,537 | `1726fedaf312889c73737ab9f2ea507927b52d28800e403e790855175aaf2c82` |

## Remaining Boundary

The prior ten-stage [native assembly integration](PUBLICATION_RUNTIME_ASSEMBLY_20261002.md)
is retained at its original controller revision. New controller composition is
unit-tested, but the current executed route is preflight only. Neither receipt
re-admits native outputs, establishes whole-base/OS compatibility, clears rights
or security, supplies a quiet window, or authorizes the 27 production timings.
The new helper/flags postdate the preserved handoff archive. Do not relabel
that correctly dated archive or rebuild it merely for this update. Scientific
source, defaults, inputs and scores remain unchanged. No DGX, contention poll,
unrelated job/service change, package upgrade or release/deposition occurs.
