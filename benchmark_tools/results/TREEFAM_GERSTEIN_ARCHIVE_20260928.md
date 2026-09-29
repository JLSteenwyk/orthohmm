# Historical Derived TreeFam Archive

The public [Gerstein lab resource page](https://info.gersteinlab.org/Ortholog_Resources)
describes an April 2009 download from TreeFam release 7 and links derived
human, fly and worm tables. Its linked archive directory remains accessible
over HTTP. HTTPS retrieval failed with an expired certificate; TLS validation
was not disabled. HTTP bytes are unauthenticated and are not admitted as
verified benchmark source material.

Downloaded the directory listing, both family-membership tables and three
gene tables (about 5.8 MB total). The [download receipt](treefam_gerstein_download_20260928.json)
records URLs, final URLs, retrieval times, response headers, sizes and SHA256
digests. Raw files remain local under `benchmarks/work/treefam_gerstein_20260928`.
The directory lists 18 derived data files, not an original NHX collection,
full MySQL dump or `treefam2reference.txt`. No downloaded code was executed.

The [content inspection](treefam_gerstein_inspection_20260928.json) found:

| File label | Parsed content |
|---|---|
| Human genes | 46,810 records |
| Fly genes | 19,789 records |
| Worm genes | 19,789 records, all data rows identical to the fly file |
| Raw family table | 8,435 families; 17,833 unique IDs; 53 absent from the downloaded gene tables |
| Multi-species family table | 934 families; 2,807 unique IDs; 52 absent from the downloaded gene tables |

The worm header says 20,151 records, but the actual file repeats the fly data
after its comments. The two family counts are each one lower than the resource
page reports (8,436 and 935). These discrepancies are observations about the
retrieved derived files, not proof of an error in the original TreeFam release.
The source page also reports historical ID-mapping problems in its own analysis.

All download hashes were checked before and after parsing; three focused
tests pass. No family or species labels from this archive were used in QfO,
and no scores, uncertainty estimates or manuscript claims changed. Even a
corrected three-species membership export would not supply the missing full
trees, their duplication events, or the benchmark-specific mapping. This lead
does not close original-source recovery. No maintainer was contacted.
