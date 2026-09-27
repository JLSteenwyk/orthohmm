# Fresh Relocated Recovery Runtime

**Subsequent portability correction:** the integrated-workflow preflight found
two unused MAFFT convenience links still pointing at the original build prefix.
Earlier native execution traces remain valid for the paths actually exercised,
but do not establish portability of those links or the entire old asset tree.
See the [fresh relative-link copy and integration test](INTEGRATED_WORKFLOW_20260927.md).
The original assets, archive and receipt have not been rewritten.

The [execution receipt](fresh_relocated_runtime_20260927.json) records a fresh
offline installation and native fixture execution outside the checkout.
Unlike the earlier harness-only relocation, this run used a new virtual
environment, copied wheels, relocated MAFFT (including helpers), relocated
FastTree, copied inputs and exported harness sources.

## Installation And Execution

`install_recovery_environment.py` refuses an existing destination, checks the
frozen lock identity and eleven-wheel inventory, creates a venv without bundled
pip, then installs with a trusted pip controller using `--no-index`,
`--require-hashes`, `--only-binary=:all:` and `--no-cache-dir`. Commands, input
hashes, logs and failures are retained; it never retries. Successful installation
commands alone do not assert a scientific or package-byte audit.

```bash
python -m benchmark_tools.install_recovery_environment \
  --base-python /absolute/python3.10 \
  --installer-python /absolute/trusted-pip-environment/bin/python \
  --wheels /absolute/copied/wheels \
  --lock /absolute/export/benchmark_tools/results/publication_recovery_requirements_20260926.txt \
  --output /absolute/fresh/environment
```

The actual new environment was under `/tmp/orthohmm-fresh-runtime-20260927`.
Its distribution inventory matches the pip report: eleven packages. All 3,052
audited installed payload files match wheel bytes before and after inference.
The existing auditor's exclusions for generated metadata/bytecode and relocated
non-site payload remain in effect; this is not exhaustive filesystem attestation.
All 69 copied wheel/tool file records match their original bytes.

The native command used isolated Python, a clean recorded environment, two
CPUs, and `MAFFT_BINARIES` pointing to the relocated helper directory. A
`strace -f -s 4096 -e trace=%file` recording followed native descendants.
Literal-prefix inspection found no reference to the original checkout/software
mount or the preceding relocated harness directory. This is bounded trace
evidence, not a filesystem sandbox or universal access proof. The base Python
installation and operating-system libraries are still shared.

## Scientific Readback

All four input files match the retained 16-gene fixture. Native execution
completed once with exit zero, and all four independent scientific readers
passed. The three root groups are byte-identical to the earlier fixture, with
36 ortholog pairs, four duplications, three speciations, two marker families,
one reconciled family and two bypassed families. No tree checkpoint was reused.
The fixture's satellite merge trace is empty. Nineteen focused installation,
export and payload-audit tests pass, including failure/timeout preservation
without retry.

## Remaining Scope

The local 99,065,669-byte archive preserves the exported harness, wheels,
relocated tools, fixture, execution trace and audit receipts; it excludes the
generated venv and temporary home directory. It is hash-recorded, not deposited
or cleared for redistribution. Independent readers ran from the original
checkout after native execution. This same-host fixture does not establish
full-data relocation, a different operating system/architecture, a rebuilt
external toolchain or publication readiness. Dedicated-host timing and the
remaining scientific and release requirements remain open.
