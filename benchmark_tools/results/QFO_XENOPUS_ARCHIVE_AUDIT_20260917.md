# Xenopus Archive Identity And Corrected Release Follow-Up

## Local Archive Verification

Read the original local2,626,624,913-byte archive without extraction. Both
selected retained FASTAs exactly match their archive members by size and
SHA-256. Iteration visited507tar members, rejected duplicate/nonregular
selected entries, and consumed gzip through EOF so its checksum/trailer was
checked. Archive and retained file hashes were checked again afterward.
No corruption or trailing-garbage exception occurred.

Archive SHA-256:
`483902f6aa9531be44ae8cb8b30eb1cb892b1c8d866ad2a1b4f96dc92e0e1978`.
Audit report `qfo_xenopus_archive_audit_20260917.json` SHA-256:
`a946c70aa8ef5cedeb58b73dc69074a3177d8e4c1485a7ac28a86e8184c97ed7`.
The independently checked canonical-to-staged byte identity remains valid.
This excludes extraction/staging loss as the explanation for the selected
files' contents. It does not authenticate the local archive against a
publisher checksum or establish that this is the intended scorer release.

```sh
python benchmark_tools/audit_qfo_source_archive.py \
  --archive qfo_benchmark/proteomes/QfO_release_2020_04.tar.gz \
  --canonical qfo_benchmark/proteomes/extracted/Eukaryota/UP000008143_8364.fasta \
  --additional qfo_benchmark/proteomes/extracted/Eukaryota/UP000008143_8364_additional.fasta \
  --output benchmark_tools/results/qfo_xenopus_archive_audit_20260917.json
```

Seven tests cover matching without extraction, missing/duplicate members,
symlinks, incorrect hashes, corrupt gzip checksum and non-gzip trailing bytes.

## Corrected Release Discovered

On2026-09-17local time, checked EBI's
[corrected2020directory](https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/previous_releases/qfo_release-2020_04_with_updated_UP000008143/)
and its [NOTES](https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/previous_releases/qfo_release-2020_04_with_updated_UP000008143/NOTES)
using HTTPS curl after the web browser tool could not retrieve the directory.
The publisher notes specify that UP000008143 alone was replaced using
2020_06 data; the other proteomes retain original2020_04 data. The directory
lists `QfO_release_2020_04_with_updated_UP000008143.tar.gz`, NOTES and README.
The corrected archive has not yet been downloaded or sequence-validated here.

The [QfO service documentation](https://orthology.benchmarkservice.org/proxy/doc)
requires a common reference-proteome set for comparable evaluation and warns
against substituting individual current UniProt proteomes because identifiers
and sequences may differ. Thus a verified, version-matched correction, not
ad hoc accession replacement, is the appropriate next investigation.

## Required Next Actions

1. Acquire the publisher's corrected archive into a new, separate location;
   retain the original archive, inputs, runs and hashes unchanged.
2. Compare corrected canonical identities and sequences with the frozen
   mapping/reference resources and the original inputs. Verify the advertised
   one-proteome scope rather than assume it from the filename.
3. Establish whether the corrected release is the sequence set used to build
   the retained2020scorer resources. Audit historical comparator input parity.
4. Prespecify any needed corrected-input reruns and label original-input
   results separately. Do not mix corrected and original inputs within a
   comparator table or silently substitute a subset-only sensitivity score.

This discovery elevates input-release matching to a publication correctness
issue. It does not prove corrected scores, change rankings, invalidate every
existing result or justify restarting active jobs on different inputs.
Current factorial work remains a frozen-original-input experiment. Preserve
it, but do not claim resolved official QfO comparability until the release
and mapping investigation is complete. Dedicated DGX scaling is unchanged.
