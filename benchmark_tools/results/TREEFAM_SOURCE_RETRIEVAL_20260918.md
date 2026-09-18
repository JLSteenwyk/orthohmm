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

Request the QfO 2020/2020.2 reference-generation source bundle from the QfO
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
