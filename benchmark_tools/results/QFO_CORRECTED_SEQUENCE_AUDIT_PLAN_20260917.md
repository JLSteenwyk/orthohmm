# Corrected Archive Sequence Compatibility Audit Queued

Job21689 runs after successful21688(file/mapping comparison), which itself
depends on21687(corrected archive acquisition). All use bizon, not the DGX.
The sequence audit reserves2CPUs/16GiB/2hours. Its committed batch script
`qfo_corrected_archive_sequences_20260917.sh` separately verifies both upstream
terminal scheduler records, final archive byte size and pinned executor state.

Executor `benchmarks/work/publication_qfo_archive_sequences_v1` is fixed at
`25394c3e2b5b7edf92b710d5edb6859b1d09866c`. The script logs Python and
Biopython versions, checks tracked analysis source cleanliness and writes to
a fresh report path. Later edits in the main worktree cannot change this job.

## What It Measures

`audit_qfo_archive_sequences.py` reads canonical sequences directly from
the archive, without extraction. It reuses the exact78-file inventory checks,
mapping coverage checks and independent native database record reader.
Every mapped input accession must agree with its numeric database entry.

Report `qfo_corrected_archive_sequences_20260917.json` separates:

- Exact sequence-byte equality.
- Byte differences completely explained by B/O/U/Z-to-X replacement.
- Unexplained sequence differences.
- Missing reference numeric IDs and unmapped input accessions.

The representation rule is explicitly restricted to the four replacements
observed in the original sequence forensic audit. It does not turn normalized
matches into exact matches or hide substitutions between standard residues.
The in-memory comparison does not change sequences on disk. Sequence agreement
among mapped entries is a distinct field from complete mapping coverage.

## Admission Boundaries

Twenty-three focused tests pass, including an archive-to-native fixture with
exact, representation-only, other sequence differences and a missing database
entry. Shell syntax passes. These establish implementation behavior, not the
corrected archive's empirical compatibility; no result is available yet.

After completion, verify scheduler status, compare archive/source hashes
between21688and21689, inspect all changed canonical files and assess both
mapping completeness and remaining sequence differences. Neither scheduler
success nor a representation match proves correct external annotations,
historical comparator input parity or improved orthology accuracy.

No corrected input freeze, extraction or inference is authorized automatically
by this job. Preserve original scores and active factorial runs; decide any
corrected-input reruns through an explicit, separately frozen protocol after
reviewing the actual compatibility evidence.
