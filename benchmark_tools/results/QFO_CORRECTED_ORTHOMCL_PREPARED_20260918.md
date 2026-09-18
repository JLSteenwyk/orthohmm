# Corrected OrthoMCL Input Preparation

Prepared native combined FASTA and genome map in the new directory
`benchmarks/results/qfo_corrected_orthomcl_v1/work`. No BLAST database,
search, inference or score has been generated. Original outputs remain
untouched. The preparation manifest is
`qfo_corrected_orthomcl_prepared_20260918.json`, SHA-256
`3106012dc42c053d42e7f9d8d08532168d5826aa6812bab4b4c232f900e3a8ff`.

An independent parser checked **984,137 of 984,137 sequences exactly**
against all 78 corrected canonical FASTAs: zero sequence differences,
complete unique IDs and correct complete species assignments in the
genome map. No residue normalization was used. This does not verify a
BLAST database that has not yet been built.

| Prepared artifact | Bytes | SHA-256 |
| --- | ---: | --- |
| all.fa | 472930503 | 11c03d6575e22c2f8bb718e59b7a06911637438c2e45181e86548647d49199c0 |
| all.gg | 23877424 | 0a3d84171661e7cb0f9cf3ebf6e739d78e80965b0e26230b319329639374d26e |

The prepared search command retains legacy BLAST 2.2.13 blastp,
E-value 1e-5, tabular format 8, 1,000 descriptions/alignments, and default
masking. It uses the full corrected combined database, with no old-release
search reuse. The proposal retains the original resource request of 180
threads, 64 later pair workers, 900 GiB and a 14-day limit; actual scheduler
allocation and environment are not yet frozen or authorized. Shared-host
resources would be descriptive, not matched timing evidence.

## Output Semantics and Failure Rules

The historical `run_orthomcl_qfo.slurm` terminates in matrix-edge conversion.
That branch must not be used for the publication comparator: the retained
publication row uses final native OrthoMCL groups expanded to cross-species
pairs, as established by the final-group scoring audit. The old launcher
also supports implicit stage reuse; the corrected run requires a separate
guarded launcher and provenance-based stage admission.

Keep legacy sequence-specific search errors visible and audit their
corrected-reference exposure. Do not disable masking, substitute another
engine, remove failed proteins from the evaluation universe, or transfer
the original 53-failure count to this new run. New failures and downstream
effects must be measured, not assumed.

## Validation and Next Steps

Twenty-three focused preparation/parity tests passed. Real preparation and
independent sequence/genome-map comparison succeeded, with source/input
hashes rechecked afterward. Next freeze configured native modules,
runtime/dependencies and a guarded formatdb/BLAST launcher, then validate
database/search/BPO/index artifacts and final native groups before scoring.
The current preparation explicitly leaves execution unauthorized.

SonicParanoid job 21710 separately passed scheduled preflight and entered
native inference. Its log confirms 78 proteomes, DIAMOND very-sensitive,
32 threads, bitscore 40, length/merging thresholds 0.75, MCL inflation 1.50
and graph-only false. Startup is not completion or accuracy admission.
