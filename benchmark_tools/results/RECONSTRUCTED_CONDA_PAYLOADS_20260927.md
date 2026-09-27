# Reconstructed Conda Payload Audit

The [audit summary](reconstructed_conda_payloads_20260927.json) checks all
6,664 declared installed entries across the 19 reconstructed base packages.
Each package archive was rehashed against acquisition pins before streaming
its payload; duplicate or unsafe member names and unexplained extra members
were rejected. Installed metadata file inventories must match the declared
paths exactly. Raw archive file digests also match the corresponding original
digests in installed metadata. No environment was modified by these checks.

| Check | Entries | Result |
| --- | ---: | --- |
| Unmodified regular-file bytes versus archive | 5,182 | Match |
| Symlink target strings versus archive | 1,286 | Match |
| Prefix-rewritten bytes versus installer-recorded digest | 195 | Match, weaker evidence |
| Regular-file bytes differing from archive | 1 | Preserved Python bytecode-cache difference |

The 195 rewritten files have not had their prefix transformations independently
reconstructed. Their installed digests agree with the installer records,
which is not equivalent to direct byte equality with the original package.
All regular-file identities and referenced metadata were rechecked afterward.
The audit covers declared Conda entries, not every extra generated file or
the separately audited pip overlay.

## Bytecode Difference

The sole direct byte mismatch is Python's
`lib/python3.10/__pycache__/_sysconfigdata__linux_x86_64-linux-gnu.cpython-310.pyc`.
It is 41,183 bytes installed versus the archive's declared 58,871 bytes.
The raw mismatch is retained in the full audit; it was not silently excluded.

Using the reconstructed interpreter, compiled the installed `.py` source into
a separate output file outside the runtime, without importing/executing that
source or overwriting the cache. The resulting file has matching 16-byte
headers and all 16 compared code-object fields, including instructions,
constants, names, filename and line table. Code-object equality also passes.
Serialized bytes still differ at 28 positions. Thus the retained cache is
consistent with compilation of the installed source under these checks, not
byte-identical to either the archive cache or the fresh serialization. This
does not prove precisely when or how the existing cache was regenerated.

## Archive-Only Notices

An initial audit assumption rejected extra payload members before creating
an output. Inspection identified 17 regular `info/licenses/` entries that
are not installed into the prefix. The completed audit inventories their
package, member path, size and SHA256 separately rather than treating them as
missing runtime files. No other extra non-directory payload was allowed.
This is not a complete notice export, component attribution or legal review.

The full audit remains at
`/tmp/orthohmm-base-reconstruction-20260927/conda_payload_audit.json`;
the summary pins it and both cache-comparison receipts. The existing native
fixture remains unchanged, and full job 22337 was not accessed or restarted
by the payload audit. Remaining transformation checks, runtime closure,
full-dataset admission, controlled timing and publication gates stay open.
