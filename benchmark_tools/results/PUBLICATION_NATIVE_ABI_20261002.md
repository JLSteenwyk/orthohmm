# Declared Native ABI Requirements

Status: 2 October 2026. Executed static inspection of the exact relocated
runtime assembly, not a scientific rerun, installation, cross-host test,
license clearance or controlled timing admission.

## Finding

All 32 compiled MAFFT helpers require the version name `GLIBC_2.34` from
`libc.so.6`. Merely satisfying the wheels' manylinux tags or supplying Linux
x86-64/AVX2 does not describe the complete assembled tool requirements.
The loader/provider selected on another host must satisfy those names;
a maximum version alone does not prove that it will do so.

The supplied FastTree has no program interpreter or GNU version-needs entries.
This is not evidence that its statically included libraries, kernel/system
requirements or license/source obligations are absent. Two other installed
MAFFT helpers are Perl scripts and are not ELF objects:
`mafftash_premafft.pl` and `seekquencer_premafft.pl`. Their interpreter and
possible external-service requirements are outside this ABI scan.

## Exact Scope

The scanner first and last validates the unchanged 102-file/18-link assembly
against its external index SHA-256
`6a496125794dca84423a0be91ab327003e15955bda91835d5ee72f28ba2839a7`.
Its original [assembly/integration admission](PUBLICATION_RUNTIME_ASSEMBLY_20261002.md)
is reused, not repeated. Native asset bytes, locks, defaults and scientific
results are unchanged. Source commit is
`48a69a833477a6824032534b25239b733a24fb63`.

| Inspected content | Count |
| --- | ---: |
| Unique wheel bytes across inference and reader roles | 12 |
| Wheel ELF objects, detected by magic rather than filename | 60 |
| Compiled MAFFT helpers | 32 |
| Supplied FastTree binary | 1 |
| Total inspected ELF objects | 93 |

All objects declare x86-64. The only observed program-interpreter path is
`/lib64/ld-linux-x86-64.so.2`. The full per-object report retains ELF class,
endianness, OS/ABI, type, flags, interpreter, library/version requirements and
raw inspector output, including weak-version flags. It distinguishes version
definitions from requirements. No version names are collapsed into a purported
universal host minimum. The three project kernels additionally require
`GOMP_4.5`, `GOMP_4.0`, `GOMP_1.0` and `OMP_1.0` from `libgomp.so.1`.

Inspection executes GNU readelf 2.42 itself, not any code from the inspected
assets. Its binary identity is recorded. The report's
`native_code_executed: false` refers to reviewed asset code, not the inspector.
The [GNU readelf documentation](https://sourceware.org/binutils/docs/binutils/readelf.html)
describes the header, segment and version-section inspection options used.
No `ldd`, import, dynamic loading or native test invocation is used.

## Validation

The new ABI mode is opt-in for the existing wheel scanner; default dynamic-tag
results remain unchanged. Parser checks reject incomplete headers/segment
tables, missing interpreter/version evidence, malformed requirements,
count mismatches, duplicate fields/versions and inspector diagnostics. The
assembly route retains external-anchor and post-read byte/mode/link validation.

All 65 focused wheel/ABI/assembly tests pass in 1.71s, zero failures/errors/skips.
An earlier overlapping 64-case panel passed before adding the incomplete
program-table check; do not sum counts. The final 8,913-byte JUnit receipt is
`benchmarks/work/publication_native_abi_validated_tests_20261002.xml`, SHA-256
`93d9b90158802ada8466389f02a61adc628d7c89b337004fa2b24e40e07c601c`.

Separate readback checks all wheel/tool ELF membership and bytes, inspector
and committed source identities, retained stdout/structured-field agreement,
the 32 helper requirements and final JUnit. It reuses the same parser and is
not an independent ELF implementation. No readelf/native workload is repeated.

| Evidence | Bytes | SHA-256 |
| --- | ---: | --- |
| [Executed inventory](publication_native_abi_20261002.json) | 1,227,076 | `75caf90467ee756d936f34260e8f2b7faff44a19774b95c74140535921191a53` |
| [Readback](publication_native_abi_validation_20261002.json) | 3,353 | `9184d206764736e4ebf11b86645f6f61e135b0de37b4abb7cb65c5ddebc480f8` |

## Boundary

This closes the selected assembly's missing declared-ABI inventory, not actual
loader resolution, CPU instruction requirements, static/dlopen/unversioned
closure, base Python/Conda/OS inspection, security, source correspondence or
redistribution review. Historical wheel/tool inventories and dated handoffs
remain unchanged. A fresh source revision is needed to obtain the new scanner;
the older validated archive is not silently relabeled as containing it.
No host contention poll, DGX action, new environment, dataset outcome or
controlled timing is introduced. The publication goal remains incomplete.
