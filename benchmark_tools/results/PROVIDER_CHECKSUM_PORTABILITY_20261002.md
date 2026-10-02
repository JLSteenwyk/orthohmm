# Provider Checksum Portability

## Defect And Correction

The actual source-548 macOS Python 3.11 fast log fails the offline Three
Kingdoms provider audit because `/usr/bin/sum -r` returns nonzero. Apple treats
the flag as a filename, not an algorithm selector. Its
[implementation](https://raw.githubusercontent.com/apple-oss-distributions/file_cmds/main/cksum/cksum.c)
selects historic BSD algorithm 1 when invoked as `sum` without flags. The
[Apple manual](https://raw.githubusercontent.com/apple-oss-distributions/file_cmds/main/cksum/cksum.1)
specifies a 16-bit rotating checksum and rounded 1,024-byte blocks. The
[GNU manual](https://www.gnu.org/s/coreutils/manual/html_node/sum-invocation.html)
documents the same default and identifies `-r` as redundant for GNU sum.
These sources support command semantics, not proof of the exact CI binary.

Change the auditor to invoke the retained `/usr/bin/sum` with just the
absolute input path. Do not replace the native algorithm, add platform
dependencies or silently retry another implementation. Record the command and
provider expectation before execution and retain stdout even when malformed.
Require three output fields, a decimal 16-bit checksum, the exact filename
and block count consistent with retained bytes. Valid differing checksums
remain `mismatch`; malformed output becomes `unresolved`, never a match.

Nonzero command exits retain return code/stdout/stderr as an unresolved row.
Execution I/O failures are also unresolved; later rows and final input/binary
checks still run. Mutated inputs continue to fail final admission. Existing
output directories cannot be overwritten. Fresh-download size/time/hash checks
and moving-release refusal remain intact. Historical provider reports and
their old auditor/command identities are not rewritten or repinned.

## Validation

**138 local tests pass in 4.51s**, zero failures/errors/skips. Scope is 32
provider cases, 15 source cases, 15 independent pair-count cases, 49 standalone
Swiss checker cases and 27 bundler cases. The earlier 32-case pass overlaps.
Sixteen provider cases are newly added. No raw source download or biological
inference/scoring occurs; network responses are offline fixtures.

Five actual native cases cover empty input, `abc`, and 1,023/1,024/1,025-byte
block boundaries, including filenames with spaces. On tested Linux/GNU
coreutils 9.4, default and explicit `-r` output agree in all five cases. The
offline audit uses the fixed known `abc` vector 16556/1, not a provider value
manufactured from the implementation under test. Additional cases distinguish
mismatches, malformed fields, wrong filename/block size, permission/nonzero
failures, no retries and postflight mutation.

[Machine-readable pins and observations](provider_checksum_portability_20261002.json)
identify source/test/JUnit files and the prior CI logs. The actual source-00dd
macOS Python 3.10 fast log confirms all 76 prior checker/bundler cases pass
(27 bundler plus 49 checker). Overall it has 13,841 passes, 29 failures,
110 skips and 30 warnings in 320.64s, including the old `sum -r` failure.
The changed auditor is not in that source, so no new-checksum macOS confirmation
follows. At 06:35:48 UTC the 3.10/3.11/3.13 test jobs are terminal failures,
full and 3.12 are live, and wheel/docs pass. Only the 3.10 log is inspected;
do not infer sibling test contents or restart incomplete handles.

## Boundaries

The eleven fixed-release downloads already verified on 2026-09-28 are reused,
not reacquired. The moving source remains unresolved. BSD sum is weak provider
metadata, not cryptographic equality; retained SHA256 and separate fresh-download
equality keep their distinct roles. This fix does not establish historical
per-tool input consumption, upstream redistribution rights, new accuracy,
controlled timing or a complete release. Frozen scientific core and scores,
prior review/statistics archives and unrelated outputs remain unchanged.
