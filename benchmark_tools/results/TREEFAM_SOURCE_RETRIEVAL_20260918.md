# TreeFam Source Retrieval Investigation

## Outcome

The requested original TreeFam-A release-7 NHX collection and
`treefam2reference.txt` have **not been retrieved**. Public reference material
was downloaded and checked, but it is not a substitute for those source
files. No family-level TreeFam uncertainty result is enabled by this work.

Downloads are isolated in `benchmarks/work/treefam_source_search_20260918/`:

| Download | Bytes | Verification |
| --- | ---: | --- |
| `ReconciledTrees_TreeFam-A.drw` | 1,998,195 | Matches Zenodo published MD5 and retained scorer byte-for-byte |
| `qfo2016_supplementary_software.zip` | 415,707 | All ZIP CRC checks pass; 421 entries; no NHX files or `treefam2reference.txt` in inventory |
| `qfo_2020_2_zenodo.json` | 10,303 | Publisher metadata names version 2020.2 and its complete file inventory |
| `treefam_home.html` | 26,572 | Snapshot of live retirement notice; informational, not a dataset |
| `treefam_download.html` | See retained file | Live download page lists TreeFam 9, not the required release 7 |

The pooled reference SHA-256 is
`6f419f96886e5cf14ac889a1437cdb2655140cb8fd281f6bb5f8e0baa69d23b3`;
published MD5 is `87202bc69e46966012cbeea456934639`.
The software ZIP SHA-256 is
`c87c17fe35024f3d06046cf0a9921b7dda99413d98d16e4cf2c8527f797e5c6f`.
The Zenodo metadata snapshot SHA-256 is
`2231f741bf41981c0d9ffaee8cf4f79baa9926ba85970f175f512d1e06f2798d`.
No raw downloads are committed into the source repository.

## Verified Sources

- [QfO 2020.2 deposit](https://zenodo.org/records/15087752),
  [machine-readable metadata](https://zenodo.org/api/records/15087752),
  [pooled TreeFam reference download](https://zenodo.org/api/records/15087752/files/ReconciledTrees_TreeFam-A.drw/content).
  Its 25-file inventory has the pooled reference, not the original NHX
  collection or `treefam2reference.txt`. Metadata specifies CC-BY-4.0.
- [2016 QfO paper](https://www.nature.com/articles/nmeth.3830) identifies
  TreeFam-A version 7 as the reference-tree source. Its
  [supplementary software](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.3830/MediaObjects/41592_2016_BFnmeth3830_MOESM332_ESM.zip)
  supplies code but not the missing source data. This older paper does not by
  itself authenticate the exact 2020 mapping file.
- The retained `generateData/AddReconciledTree.drw` explicitly reads
  `../data/treefam/treefam2reference.txt`, iterates `../data/treefam/*.nhx`,
  and unions their relations into one `TreeFamA` case. The currently exposed
  GitHub tree does not contain those source assets.
- [TreeFam download page](https://www.treefam.org/download) exposes release-9
  data, which must not silently replace release 7. The
  [live homepage](https://www.treefam.org/) announces service retirement on
  September 30, 2026 and directs archival questions to
  [EMBL-EBI support](https://www.ebi.ac.uk/about/contact/support), topic
  `TreeFam - database of animal gene trees`.

## Unsuccessful Retrievals

These are observations from this search, not proof that no copy exists:

- HTTPS `ftp.sanger.ac.uk/pub/treefam/` and `pub2/treefam/`: HTTP 404.
- FTP `ftp.sanger.ac.uk/pub/treefam/release-7.0/`: connection timeout;
  not evidence that the archive itself is absent.
- HTTPS EBI `pub/databases/treefam/` and `pub/databases/TreeFam/`: HTTP 404.
- HTTPS TreeFam `static/download/release-7.0/`: HTTP 404.
- HTTPS `legacy.treefam.org`: certificate hostname mismatch; verification
  was not disabled. HTTP legacy site redirects to the release-9 homepage.
- QfO `refsets/data/treefam/` and its `treefam2reference.txt`: HTTP 404.
  The `refsets/` index returned 403; no access restrictions were bypassed.
- Internet Archive CDX query for the Sanger release-7 gzip URL pattern
  returned an empty list. A second mapping-name query timed out. Neither
  establishes global absence from web archives.
- Exact filename searches and inspected publisher/GitHub inventories did
  not reveal a downloadable original mapping.

## Required Next Step

### September 20 PhyloMCL Archive Check

A fresh public search identified the [PhyloMCL materials archive](https://sourceforge.net/projects/phylomcl/files/Materials/),
including `QfO_data.tar.gz`. Downloaded the complete archive from
`https://downloads.sourceforge.net/project/phylomcl/Materials/QfO_data.tar.gz`
to `benchmarks/work/treefam_source_search_20260918/phylomcl_QfO_data.tar.gz`.
Its SHA-256 is
`fb3a0a2b6c304fe6fd6b5221a9013c046482f4302c382d6db98c668556e4524a`.
This is a locally computed identity, not an independently published checksum.

`tar -tzf` completed successfully. The complete inventory contains FASTA
files, `gene.length`, `gene.idmap`, `66_species.nwk`, and PhyloMCL ortholog
pairs/groups, but no NHX collection, nested source archive, or
`treefam2reference.txt`. No archive contents were extracted or executed.
The archive is retained outside Git and is not admitted as a QfO 2020
reference-generation source or a substitute mapping.

The existing pooled reference and supplementary software ZIP were rehashed
and still match their recorded SHA-256 values. Fresh exact-filename searches
did not locate the missing mapping. The TreeFam homepage still directs
archive questions to EMBL-EBI support. No maintainer contact was sent; the
original trees and mapping remain missing, and benchmark scores and
family-level uncertainty claims remain unchanged.

### Additional September 19 Source Check

Re-downloaded the QfO deposit metadata to
`benchmarks/work/treefam_source_search_20260918/qfo_2020_2_zenodo_recheck_20260919.json`
(10,303 bytes; SHA-256
`6eee9d887cbb5508e703e21420cd26964b5bb4284ae3eec2b0efa9f3d7e8ce33`).
Inspected all 25 file entries: the pooled TreeFam reference is present,
but the original NHX collection and `treefam2reference.txt` are not listed.
The general `mapping.json.gz` must not be treated as the missing
TreeFam-specific reference-generation mapping.

Fresh exact-filename and release-7 searches did not identify a download.
The original [TreeFam publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC1347480/)
also identifies a historical Chinese mirror at `treefam.genomics.org.cn`.
An HTTPS request to that host timed out at connection establishment after
10 seconds; its contents could not be checked. This does not establish
that the source archive is absent.

Rehashed the retained pooled reference and supplementary-software ZIP;
both still match the checksums above. No original tree or mapping was
recovered, no substitute release was admitted, and no maintainer contact
was sent. The specific maintainer request below remains the next step.

### September 19 Retrieval Check

Fresh exact-name and release-7 searches did not recover the original mapping
or NHX collection. The [Gerstein laboratory resource inventory](https://info.gersteinlab.org/Ortholog_Resources)
documents an April 2009 download of release-7 MySQL tables and links derived
human/fly/worm tables, not the QfO mapping or original tree collection. A
request to its linked archive directory,
`https://archive.gersteinlab.org/proj/orthologs/Orthologs/TreeFam/`, failed
TLS verification because the certificate had expired. Verification was not
disabled; the directory contents remain unverified. A fresh HTTPS request
to the Sanger release-7 directory returned HTTP 404.

Both previously downloaded reference/software files were rehashed and still
match the SHA-256 values above. No additional source trees or mapping were
downloaded, and no benchmark results or uncertainty estimates were changed.

### Additional Public Leads Checked

- Repeated exact searches for `treefam2reference.txt` and
  `treefam2reference` did not identify a downloadable original mapping.
- The [TreeSoft SourceForge inventory](https://sourceforge.net/projects/treesoft/files/)
  lists TreeFam Perl API releases and tree software. Its
  [OldFiles directory](https://sourceforge.net/projects/treesoft/files/OldFiles/)
  lists three software archives, not a TreeFam release-7 dataset. These
  visible inventories do not establish what every software archive contains.
- An Internet Archive availability request for
  `ftp.sanger.ac.uk/pub/treefam/release-7.0/` with timestamp `20120101`
  returned HTTP 429. This is an unsuccessful lookup, not evidence of
  absence from the archive; no rate-limit bypass was attempted.
- Rehashed the retained pooled reference and supplementary software ZIP;
  both still match the SHA-256 values recorded above.

No additional original trees or mapping were downloaded from these leads.

### Follow-up Public Archive Search

A further exact-filename and release-7 search did not locate a downloadable
mapping or the original NHX collection. Additional leads checked:

- Downloaded the [EBI archived-database inventory](https://ftp.ebi.ac.uk/pub/databases/archived_databases_200623.txt)
  to `benchmarks/work/treefam_source_search_20260918/ebi_archived_databases_200623.txt`.
  The 41,860-byte file has SHA-256
  `a42dc78361cb8298c70af30e202d063a355b020789f05b02d17eea5711f7e982`.
  A case-insensitive search for `treefam` or `tree.fam` found no match. This
  inventory check does not establish absence from all EBI archives.
- The [Naturalis TreeFam data-mining tutorial](https://naturalis.github.io/mebioda/doc/week1/w1d5/lecture1.html)
  leads to [a public download script](https://github.com/rvosa/bh15-fossil-paralogy/blob/master/pipeline.sh).
  Its source is the same `static/download/treefam_family_data.tar.gz` exposed
  by the release-9 download page, not an independently identified release-7
  archive. It was not admitted as a replacement.
- Recomputed both retained pooled-reference and supplementary-ZIP SHA-256
  checksums; they still match the values above.

No original source file was recovered in this follow-up, no benchmark score
was changed, and no contact request was sent.

### Additional Repository-History Check

A bare clone of `https://github.com/qfo/benchmark-webservice.git` is retained
at `benchmarks/work/treefam_source_search_20260918/qfo-history.git`, with HEAD
`c0854a96c1a0fd7f2a891d971af0863002fabc90`. A filename-history search across
all fetched refs for `*treefam*`, `*TreeFam*` and `*.nhx` found the TreeFam-A
benchmark metadata and several SwissTree/example NHX files, but not the
requested TreeFam collection or mapping. This searches reachable fetched
history, not deleted upstream refs or private archives. The GitHub commits
API query for `data/treefam` also returned an empty list.

The pooled reference and supplementary-software ZIP were rehashed and still
match the SHA-256 values above. No newly recovered original TreeFam input
has been admitted for analysis. The [Sanger archive page](https://www.sanger.ac.uk/tool/treefam/)
explicitly states that the resource is no longer available at Sanger.

### OrthoFinder2 Supplement Inspected (2026-09-28)

Job 22344 completed 0:0 in 23m43s. All 1,957,180,078 bytes were retrieved;
MD5 matches the public record and independent SHA256 reproduction gives
`44d0b6825abead63ef9d7f9b956abe1d6c9cb7e90cc23e212d1ad818f54f35e4`.
The tar inventory contains 548,923 members. No member name matched TreeFam,
`treefam2reference`, `.nhx` or reference; only README and its backup matched
the candidate-name filter. No nested tar/zip/tgz filenames were identified.

Read the exact 2,456-byte README through the archive API without extracting
paths or executing archive contents. It describes QfO proteomes and uploaded
OrthoFinder/OrthoMCL prediction files, plus fungal/chordate analyses and gene
duplication simulations. The six QfO result filenames agree with that
description. This does not identify the original TreeFam-A release-7 trees
or QfO mapping. Do not repurpose inferred/simulated trees as reference family
units. This is a bounded negative finding from filenames and documentation,
not a claim that every tree or prediction file was semantically inspected.

The [retrieval receipt](orthofinder2_treefam_search_22344.json) records the
archive identity, inventory/README hashes, relevant paths and limitations.
The large archive and full inventory remain local; no raw dataset was added
to Git. The original-family uncertainty limitation remains unresolved.

### OrthoFinder2 Supplement Retrieval Started (2026-09-28)

Revisited the remaining public [Zenodo record 1481147](https://zenodo.org/records/1481147).
Its linked `OrthoFinder2_Zenodo.tar.gz` is available over verified HTTPS;
the server reports 1,957,180,078 bytes and MD5
`5907cdf22b211f2c4e6bead42eb17878`, matching the record's published checksum.
This establishes an available archive, not the presence of original TreeFam
trees or `treefam2reference.txt`.

The initial foreground transfer was deliberately stopped (exit 143) after
89,489,408 bytes to move this approximately 35-minute download to tracked
Slurm job **22344**, one CPU, 2 GiB RAM, two hours, no requeue. The scheduled
transfer explicitly resumes those bytes using HTTP range support; no second
concurrent transfer was launched. Its live state was verified as RUNNING.
The local directory is `benchmarks/work/treefam_orthofinder2_public_20260928`.
No automatic retry is configured. Inspect the same job before any resume.

`inspect_orthofinder2_archive.py` will verify exact size and published MD5,
record SHA256 and inventory tar member names without extracting or executing
anything. Candidate names remain leads; nested archives and original-release
identity still require inspection. Four unit tests passed (0.14 seconds),
including checksum/size rejection and non-extraction of symlink entries.
Originals are not yet recovered. No maintainer was contacted or DGX accessed.

### Release-7 API Archive Inspected (2026-09-27)

Downloaded the [TreeSoft release-7 Perl API archive](https://downloads.sourceforge.net/project/treesoft/TreeFam-Perl-API/7v1/Treefam-7v1.tar.gz)
to `benchmarks/work/treefam_source_search_20260918/Treefam-7v1.tar.gz`.
The 40,113-byte archive has SHA-256
`d862953bf2b4968efc538d700bec3118c5fbc0f608f163020fbafc15fbdfeec6`.
Its complete tar inventory contains nine Perl modules and Subversion
metadata/base copies, not NHX trees, a database dump, a nested data archive,
or `treefam2reference.txt`. No downloaded code was executed or installed.

Reading `Treefam/Config.pm` directly from the archive identifies API version
7 and the public anonymous database `treefam_7` on
`vegasrv.sanger.ac.uk:3308`. A local `getent hosts vegasrv.sanger.ac.uk`
lookup returned no address (exit 2); no database connection was attempted.
This unsuccessful lookup does not establish permanent database loss.
The archive therefore closes one previously uninspected software-package
lead but does not recover the benchmark inputs. No scores or uncertainty
estimates changed, and no maintainer contact was sent.

### Public-Only Search Resumed (2026-09-27)

The user explicitly requested internet searches instead of contacting anyone.
The contact draft below is historical and **must not be sent** under the current
instruction. Renewed exact-filename and release-7 searches did not recover
`treefam2reference.txt` or a verified original NHX collection.

The public [CyVerse TreeFam integration notes](https://cyverse.atlassian.net/wiki/spaces/iptol/pages/242170711/TreeFam)
document the historical anonymous `db.treefam.org:3308` database and list
`treefam_7`. They also describe tree storage, but do not supply the requested
archive or QfO mapping. Local DNS lookups for both that host and the earlier
`vegasrv.sanger.ac.uk` returned no address; no database connection was attempted.
This does not establish permanent loss. [InterMine's own documentation](https://app.readthedocs.org/projects/intermine/downloads/pdf/latest/)
points back to the already investigated Sanger release-7 MySQL path and two
gene/ortholog tables, not an independent tree/mapping mirror.

The [Selectome paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC2686563/)
identifies releases based on TreeFam 4 and 6, not the required release 7;
those derived trees cannot silently replace the benchmark source. The public
[OrthoFinder2 supplement](https://zenodo.org/records/1481147) exposes a 2.0-GB
archive, but its listing does not establish that the missing originals are
inside. It was not downloaded or declared a recovered source in this check.
The Broccoli Zenodo landing page returned HTTP 429; no rate-limit bypass was
attempted. No original source file was downloaded, admitted or substituted,
and no email, issue or support request was sent.

### Broccoli Supplement Available Again (2026-09-28)

The ordinary public [record](https://zenodo.org/records/3710751) and
[metadata API](https://zenodo.org/api/records/3710751) now respond successfully;
no rate-limit bypass was used. The metadata identifies `data_Zenodo.zip`,
719,576,680 bytes, MD5 `fd77f9ef7a5b82a87143b603902901d4`, and describes QfO
2018 results among its contents. This does not prove the presence of the
original TreeFam trees or the mapping used by the retained QfO 2020 reference.

After an initial 120-second transfer timed out with 202,557,640 bytes retained,
tracked local job 22357 explicitly resumed the partial file. Its workflow
checks complete size and MD5, computes SHA256 and inventories ZIP members
without extraction or execution. Job state was verified RUNNING. See
[submission evidence](broccoli_treefam_submission_22357.json). Content
inspection and checksum validation remain pending; no originals are admitted.

### Historical Contact Draft (Not Authorized)

The earlier proposed next step was to request the QfO 2020/2020.2 reference-generation source bundle from the QfO
maintainers, and the release-7 archive from TreeFam/EBI if necessary. A draft
request follows; **no email, issue or support submission has been sent**.

> Subject: Original TreeFam-A inputs for QfO reference dataset 2020.2
>
> We are auditing the TreeFam-A benchmark in the QfO 2020.2 reference deposit
> (Zenodo 15087752). Could you provide the individual TreeFam-A NHX trees and
> `treefam2reference.txt` used by `generateData/AddReconciledTree.drw`, along
> with the release/build identity, mapping-generation procedure, checksums
> and reuse license? A persistent source archive would be ideal.
>
> Our pooled `ReconciledTrees_TreeFam-A.drw` matches the deposit MD5
> `87202bc69e46966012cbeea456934639` (1,998,195 bytes). We need the original
> family identities and mapping to validate family-level resampling, not a
> newer TreeFam release or a new similarity-based mapping. Please include
> the accepted/rejected family list and relation-extraction logs if retained.

Before using recovered files, reproduce the pooled membership and labeled
relation union, check overlapping families/conflicting labels and isolated
members, and establish a defensible resampling unit. Do not equate connected
components of the pooled relation graph with original source families.
