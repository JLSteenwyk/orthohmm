# Xenopus Provider Follow-Up

The [bounded public-provider receipt](xenopus_provider_metadata_20260930.json)
retains the successful official archive-directory and README responses, the
individual historical-file 404, and current-file HEAD metadata. Four local
compressed/raw/staged/BUSCO file pins were checked before and afterward.
No retained input or benchmark score changed.

The [official 2025_04 knowledgebase archive listing](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2025_04/knowledgebase/)
offers a 140 GB knowledgebase archive, a documentation archive and a Swiss-Prot
archive, not an individually listed reference-proteome directory. The
[archive README](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/README)
describes the knowledgebase archive as DAT files and documentation; it is not
a promise of reproducing the historical per-proteome FASTA serialization.
The constructed individual 2025_04 proteome URL returned HTTP 404.

The retained gzip is 14,143,121 bytes. Current HEAD metadata reports 14,105,984
bytes and Last-Modified 3 September 2026 for the moving source URL. This is
metadata evidence of different compressed size, not a full raw-sequence
comparison. The retained gzip header's MTIME is 25 November 2025; it is not
proof of a UniProt release or download date. No historical release is assigned.

This closes neither the historical-release question nor every possible public
mirror. No 140 GB download, fresh proteome replacement or inference rerun was
attempted. The earlier eleven fixed Ensembl reacquisitions remain valid and
unchanged. Preserve this concrete negative lead to avoid repeating these same
URLs on each resumption; continue only with a new candidate or changed evidence.
