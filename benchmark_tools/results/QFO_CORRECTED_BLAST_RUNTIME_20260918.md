# Corrected OrthoMCL Legacy BLAST Runtime Inventory

Recorded and reverified 67 runtime-tree entries covering the installed
BLAST 2.2.13 distribution (executables, matrix/data files and documentation)
and the resolved host libraries observed with `ldd`: libm, libpthread,
libc and the ELF loader. The loader symlink and resolved target checksum
are recorded. Inventory SHA-256:
`9ce46b5329f34384b03980b62fbfe9522f244dc227ce6506071ce7000b154e42`.

`ldd` reported no missing library for blastall or formatdb. The distribution
VERSION file records a December 2005 x86_64 build; that file alone is not
an executable version probe. No `.ncbirc` was present in the current
repository root or user home when inspected. These absence checks are
observations, not future-run guarantees.

The inventory is prospective identity evidence, not a historical execution
trace, complete operating-system snapshot, or inference authorization.
The guarded launcher must still freeze the effective environment,
configuration/data lookup, actual scheduler allocation and stage commands;
recheck absence/presence conditions and inventory before/after execution;
and validate formatted database/search artifacts separately. No formatdb
or BLAST inference was run in this milestone, and no legacy defaults or
failed-query handling were changed.
