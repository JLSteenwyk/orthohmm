# Live Search-Stage Library Snapshot

The [snapshot receipt](integrated_search_libraries_20260927.json) records
read-only observations of the native parent PID 2642246 and one search worker
PID 2642260 in job 22337. Scheduler membership was checked before and after;
both command lines exactly match the recorded native stage command. Process
start ticks and command lines remained unchanged across collection. Raw
`/proc/PID/maps`, process metadata, backing-file identities, scheduler output
and the collection command are retained. No environment variables or secrets
were collected. No debugger, signals, injection or package installation was
used, and no native workflow was rerun.

## Observed Libraries

| Scope | Mapped library paths | Exact selected-wheel member matches |
| --- | ---: | ---: |
| Parent | 67 | 23 |
| Sampled search worker | 72 | 26 |
| Union | 72 | 26 |

The sampled worker maps `hmm_viterbi.so` and `kmer_prefilter.so` whose backing
files match the frozen OrthoHMM wheel members exactly. It additionally maps
Numba's `omppool` extension and `/home/bizon/anaconda3/lib/libgomp.so.1.0.0`.
No TBB library appears in either snapshot. This is evidence for the observed
search stage, not proof that TBB or other libraries are unused in every
worker, optional path or later stage.

The 46 paths without selected-wheel byte matches belong to the base Anaconda
Python installation or `/usr/lib/x86_64-linux-gnu`. They include standard-
library extensions, compression/crypto libraries, C/C++ runtimes and the
system loader. The selected wheel lock therefore does not by itself capture
the complete runtime used here. External MAFFT/FastTree processes and later
clustering/reader stages are not covered by this snapshot.

All 139 per-process backing-file records were rehashed independently after
collection. Hashes describe files on disk, not in-memory code identity or
complete static-library attribution. Basename/SONAME candidates from the
earlier ELF inventory are not substituted for actual file-byte matches.
Package source/notice attribution and cross-host restoration remain open.

## Collection Correction And Limits

An initial collection command failed before writing its output: it appended
candidate metadata directly to the file-record dictionary, then passed that
enriched dictionary to the existing exact-record checker. The checker compares
the entire dictionary, so it rejected the extra field; this was not evidence
that a library changed. The successful observation nests the unmodified file
record separately from candidate annotations. No prior artifact was overwritten
and no inference restart occurred.

The existing `mapped_libraries` parser was reused. It preserves spaces in
paths and rejects deleted libraries; it selects names containing `.so`, not
every possible mapped executable. The Python executable is recorded separately.
This observation does not admit job 22337, demonstrate biological accuracy,
establish controlled timing, clear security/license obligations or prove
publication readiness.
