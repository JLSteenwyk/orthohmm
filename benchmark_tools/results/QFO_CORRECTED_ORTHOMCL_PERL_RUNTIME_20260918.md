# Corrected OrthoMCL Perl Runtime Inventory

## Snapshot

`qfo_corrected_orthomcl_perl_runtime_20260918.json` records 30,081 entries
(8,981,638 JSON bytes), including the complete 822-MiB native environment,
OrthoMCL script/module, MCL executable, selected system libraries, shell and
`date`. The manifest SHA-256 is
`9cc29e84777f47064c23096ce36e32d5465616979b05febe510219d58b6b3239`.
Only metadata/hashes are committed; no runtime binaries or raw datasets are
added to the repository.

The three external symlinks are `/bin/sh`,
`/lib/x86_64-linux-gnu/libcrypt.so.1`, and
`/lib64/ld-linux-x86-64.so.2`. Their resolved regular-file sizes and hashes
are captured by the inventory. There are no other external symlinks in this
snapshot. Interpreter and MCL linkage was inspected locally to select the
system-library roots; this is not a complete operating-system image.

## Observed Module-Load Probe

`inspect_orthomcl_perl_runtime.py` verifies the complete snapshot, launches a
native Perl probe in a fresh directory with the clean environment, checks
every observed module and mapped file against the inventory, then verifies
the complete snapshot and checked artifacts again.

The successful v2 probe bound 100 loaded modules and 15 mapped binary/library
files. Observed versions:

| Component | Version |
| --- | --- |
| Perl | v5.26.2 |
| BioPerl SearchIO | 1.007002 |
| Storable | 3.15 |
| OrthoMCL | 1.4 |

The report is `qfo_corrected_orthomcl_perl_probe_20260918.json`; raw JSON/log
outputs are in `benchmarks/work/orthomcl_perl_runtime_probe_v2_20260918`.
The first failed probe is preserved in the non-v2 directory.

## Current-Directory Lookup Finding

This Perl build appends `.` to `@INC`, including when directly tested with
`PERL_USE_UNSAFE_INC=0`. The first strict probe therefore failed rather than
concealing this lookup. The reviewed v2 inspector accepts this observed
default only when its fresh working directory contains no files other than
its two ordinary, nonsymlink output files. All actually loaded modules must
still match the snapshot. Other relative paths or executable lookup hooks
are rejected.

This isolated module-load probe does not approve unrestricted `.` lookup in
production. The upcoming launcher must explicitly constrain its working
directory or remove that lookup without changing native scientific behavior.
Configured native sources and pair-parallel modifications still require their
own exact source freeze and validation.

## Verification And Limits

All 17 runtime-inspector/inventory tests passed, covering changed versions,
unknown external links, unbound/changed loaded files, relative/hook paths,
empty inventories and unexpected/symlink working-directory files. The real
v2 probe passed before/after inventory checks.

Shell, date and MCL are inventoried but are not executed by this Perl probe.
The observed mappings do not establish every future import or syscall.
Before/after equality cannot exclude temporary changes during execution.
Python conversion runtime, configured scripts, job provenance and full
conversion/index outputs remain separate integration requirements.

No search or inference job was started, changed or restarted. Both snapshot
and probe reports leave scientific execution and accuracy unauthorized; no
corrected benchmark result or publication-ready claim follows from them.
