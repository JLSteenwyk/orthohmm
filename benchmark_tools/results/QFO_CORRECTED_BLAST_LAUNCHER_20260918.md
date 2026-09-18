# Corrected Legacy BLAST Launcher

The launcher binds the prepared OrthoMCL manifest
`3106012dc42c053d42e7f9d8d08532168d5826aa6812bab4b4c232f900e3a8ff`
and legacy BLAST runtime inventory
`9ce46b5329f34384b03980b62fbfe9522f244dc227ce6506071ce7000b154e42`.
It runs only fresh formatdb and blastall stages. BPO conversion, clustering,
final-group pairs and scoring remain behind separate admission gates.

The search retains the frozen blastp command, E-value 1e-5, format 8,
1,000 descriptions/alignments and default masking; no failed queries are
pre-removed. The environment is explicit and excludes inherited matrix,
NCBI and loader overrides. The runner rejects newly present user/working
directory NCBI configuration or checked global configuration/preload files.
It does not install a replacement matrix or configuration.

It verifies prepared inputs, original corrected source records, runtime
and exact commands before execution and after each stage. Existing search
artifacts are refused. Formatted database files must exist and be nonempty
before search; their hashes must remain unchanged afterward. Successful
BLAST output is retained pending full database/query/log admission, not
silently treated as complete orthology inference. All command failures
retain status, logs and partial artifacts.

Validation: real preflight passed without starting search; 19 focused
tests passed including an opt-in installed-engine smoke. The smoke built
a tiny database and recovered all four expected self/cross hits for two
identical test sequences using the explicit environment and legacy defaults.
It verifies basic runtime/configuration viability, not scientific accuracy,
large-database completeness or multithread equivalence. Shell syntax passed.

Native search may be submitted from a committed pinned executor with the
original proposed allocation: 180 CPUs, 900 GiB, 14-day limit on bizon,
no requeue. It may wait for sufficient resources; no unrelated job should
be stopped or resource gate bypassed to start it. This is descriptive
shared-host accuracy work, not part of the dedicated DGX timing series.
