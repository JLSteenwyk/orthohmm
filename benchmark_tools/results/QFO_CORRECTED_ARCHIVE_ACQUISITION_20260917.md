# Corrected QfO Archive Acquisition And Comparison Preparation

The prior archive audit established that the retained original canonical and
additional Xenopus FASTAs match the local original archive. The publisher
provides a distinct corrected archive. Acquisition is isolated under
`benchmarks/work/qfo_corrected_source_20260917/`; no frozen input is replaced.

Source: [EBI corrected archive](https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/previous_releases/qfo_release-2020_04_with_updated_UP000008143/QfO_release_2020_04_with_updated_UP000008143.tar.gz).
HTTPS HEAD on2026-09-18UTC returned200, Content-Length2,648,666,198,
ETag`"9ddf7056-5b478e2a3e88c"`, Last-Modified19November2020 17:16:15GMT
and Accept-Ranges bytes. These are server metadata, not a publisher SHA-256.

The initial foreground transfer was deliberately interrupted with SIGINT
(exit130) to move slow acquisition into a durable scheduler job. This was
not an observation timeout or restart of scientific computation. Its
79,556,608-byte partial file was retained. Slurm job21687 on bizon resumes
that same file with curl `--continue-at -` and an If-Match ETag guard.
It uses1CPU/2GiB, never the dedicated DGX. The committed batch script
`qfo_corrected_archive_acquire_20260917.sh` verifies final length and prints
SHA-256, but successful acquisition alone will not validate archive contents.
Log: `benchmarks/work/qfo_corrected_source_20260917/acquire_21687.log`.

## Prepared Comparison

`compare_qfo_corrected_archive.py` is implemented and seven tests pass.
It streams every canonical FASTA from the corrected archive, checks the
exact78-file inventory against the frozen original manifest, records changed
file hashes, counts unique sequence accessions and numeric mapping coverage,
and inventories recovery of the14missing SwissTrees accessions. It rejects
duplicate file/accession/numeric identities, unexpected or missing proteomes,
malformed headers and incomplete gzip reads. Additional/DNA files are not
canonical inputs. It performs no extraction or input replacement.

After scheduler completion and size/hash inspection, run from repository root:

```sh
python benchmark_tools/compare_qfo_corrected_archive.py \
  --archive benchmarks/work/qfo_corrected_source_20260917/QfO_release_2020_04_with_updated_UP000008143.tar.gz \
  --prepared benchmark_tools/results/qfo_factorial_prepared_20260917.json \
  --mapping qfo_benchmark/benchmark-webservice/reference_data/2020/mapping.json.gz \
  --aliases benchmark_tools/results/swiss_sequence_alias_audit_20260917.json \
  --output benchmark_tools/results/qfo_corrected_archive_comparison_20260917.json
```

No corrected-release comparison result exists at this milestone. Do not
infer missing-protein recovery, one-proteome-only change, scorer compatibility
or improved accuracy from the publisher filename or passing synthetic tests.
Those require the completed empirical audit. Preserve original-input scores
and ongoing factorial runs while this correctness investigation proceeds.

## Pinned Dependent Comparison Queued

Submitted comparison job21688 with `--dependency=afterok:21687` on bizon,
2CPUs/16GiB/2hours. Its script is
`qfo_corrected_archive_compare_20260917.sh`. The script separately requires
the exact successful acquisition scheduler record and final archive byte
count, so dependency submission itself is not evidence of a valid download.

The executor is detached worktree
`benchmarks/work/publication_qfo_corrected_archive_audit_v1`, pinned to
`a1fbffcc8cb5b4fdabca18843d79d94b0baadaef`. It checks HEAD and tracked analysis
source cleanliness before execution, records Python/Biopython versions in
its log and uses the extended all-numeric-ID/species-interval comparison.
The combined comparison/archive tests pass20cases; shell syntax passes.

Output remains `qfo_corrected_archive_comparison_20260917.json`, written
only to a fresh path after the complete audit returns. The frozen-input
manifests, archive and mapping are checked before/after reading. No input
extraction, replacement or inference is automatically scheduled. Independently
inspect terminal status, changed-file inventory, missing IDs, sequence
compatibility and input/scorer versions before declaring the release suitable.

The retained `reference_data/2020/ServerIndexed.db` contains sequence entries
as well as mapping identifiers. Investigate direct sequence-content comparison
with that resource: accession coverage alone is insufficient to establish
sequence-byte compatibility. This stronger check is not yet implemented or
claimed complete by job21688.
