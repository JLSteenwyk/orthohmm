# Legacy BLAST Recovery Source Evidence

## Recovered Release

The earlier `legacy/` URLs were obsolete. NCBI's executable index points
to [legacy.NOTSUPPORTED/2.2.13](https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.13/).
Downloaded `ncbi.tar.gz`, `blast-2.2.13-x64-linux.tar.gz`, and `MD5SUM.txt`
from that exact directory on September 23. Retained outside Git under
`benchmarks/work/legacy_blast_source_20260923/`. No downloaded executable
was installed or run; selected source files were extracted for reading.

Both archive MD5 values match NCBI's published file (historical integrity
checks, not a modern cryptographic authenticity guarantee):

| File | MD5 | SHA-256 |
| --- | --- | --- |
| `ncbi.tar.gz` | `c42e3360ea7fcbb8b17bdaa5fc8f79ed` | `3f311eb066a49c73eef36cdb8305c97c18ed35f7f9bc2f8904e18d4298e5a85d` |
| `blast-2.2.13-x64-linux.tar.gz` | `c45e01ad09cec17b9af338df53c8b07f` | `e284b4f95adf52267b21d3cbdffb9987c8ef038c557e0b3527e54fa97796bb0c` |

`MD5SUM.txt` SHA-256:
`8f6c4abdf3396fef4b0e11d58ecfb0450f2954b68296d40efc707cc2fca09948`.
Extracted `ncbi/demo/blastall.c` SHA-256:
`3344597098b9ce248be14cab885646a9a9604c7e6021f143b4ede356e172ff60`.

The archived `blast-2.2.13/bin/blastall`, read directly from the tar file,
and the installed `/SOFTWARE/blast-2.2.13/bin/blastall` beneath the project
disk both have SHA-256
`34b91fedd2e7858a3478e83758e3f5339c09f96f543d9338488ea4b84628fb75`.
This establishes installed binary identity with the official release,
not a reproducible source-to-binary build.

## Inspected Control Flow

Line numbers refer to that archived `ncbi/demo/blastall.c`:

- Line 959 defaults `ARG_FORCE_OLD` (`-V`) to true. Lines 2602-2611
  choose `Main_old` unless the new-engine conditions are satisfied.
  The frozen command did not override `-V`.
- Lines 1799-1810 obtain `NumQueries` and reject concatenation for
  programs other than blastn/tblastn. The actual run is blastp without `-B`.
- The main loop starts at 1837. Lines 1900-1940 read one FASTA entry
  per iteration for this non-megablast, non-concatenated path.
- Line 2099 calls `BioseqBlastEngineWithCallbackMult` before results are
  formatted. Lines 2213-2226 print tabular results for this query before
  proceeding through cleanup and the next iteration (loop ends at 2508).
- The source contains output flushing, but this does not establish
  filesystem durability at an abrupt interruption.

These observations support a sequential query-output hypothesis for the
retained run. They do not by themselves validate every byte in the partial
file, prove its final query complete, or authorize skipping no-hit queries.

## Recovery Work Still Required

1. Inventory query blocks and FASTA ordering across the retained table,
   validating all retained complete rows with the existing numerical and
   alignment checks. Preserve the original bytes and record the proposed
   prefix boundary separately; reject noncontiguous or out-of-order blocks.
2. Exclude the entire last query block, not just its malformed row. Treat
   every query without trustworthy retained results as requiring replay,
   including no-hit/failed queries when completion cannot be established.
3. Replay a prespecified small set of earlier queries and boundary queries
   against the unchanged complete database with the same binary, scoring,
   filters, hit limits, and environment. Compare per-query HSP results and
   query diagnostics. This is an accuracy-equivalence check, not timing.
4. If evidence supports reuse, freeze a separate recovery plan with an
   exhaustive query partition, provenance, output merge rules, and durable
   per-batch completion records. Verify full coverage and no duplicated
   query blocks before downstream admission. Otherwise run a separate
   fresh search. Do not append to or relabel the interrupted output.

No recovery search or output merge has been run. Job 21746 remains held.
This report resolves the source-location question, not recovery admission.
