# Five Base Runtime Files Verified Against Archives

The [acquisition and comparison receipt](integrated_base_archive_comparison_20260927.json)
resolves the five installed-file digest gaps in the earlier
[base-runtime attribution](INTEGRATED_BASE_RUNTIME_PACKAGES_20260927.md).
Each exact package archive was downloaded over verified HTTPS from its
recorded `repo.anaconda.com/pkgs/main/linux-64` URL. Every archive matches
the SHA256 retained in local Conda metadata. Independent `sha256sum` output
also matches those five pins. Historical metadata and receipts are unchanged.

| Package | Archive Bytes | Selected Library | Installed Bytes Match |
| --- | ---: | --- | --- |
| bzip2 1.0.8 h7b6447c_0 | 80,125 | libbz2.so.1.0.8 | Yes |
| libgcc-ng 11.2.0 h1234567_1 | 5,602,184 | libgcc_s.so.1 | Yes |
| libgomp 11.2.0 h1234567_1 | 485,145 | libgomp.so.1.0.0 | Yes |
| libuuid 1.41.5 h5eee18b_0 | 28,110 | libuuid.so.1.3.0 | Yes |
| zlib 1.2.13 h5eee18b_0 | 105,851 | libz.so.1.2.13 | Yes |

Local downloads are retained under
`benchmarks/work/publication_base_runtime_archives_20260927/`, not committed
or installed. The receipt stores exact URLs, archive and metadata identities,
selected tar member/type/size, decompressed member SHA256, installed-file
identity and the executed collection command. ZIP member uniqueness and one
package payload were required; the payload was streamed through zstandard
0.19.0 and the standard tar parser. Only the exact regular-file member was
hashed. Nothing was extracted into the environment, and no package code ran.
All input, archive, metadata and installed-file records were rechecked after
comparison. No retry, package upgrade or inference restart was needed.

## Interpretation

The preceding local inventory had 35 matching Conda installed-file digests,
seven matching system-package MD5 entries and five unknown Conda file digests.
The latter five now have direct byte equality to hash-pinned package members.
This is stronger than inferring identity from a package owner or version.
It is still limited to previously observed search-stage files and the Python
executable; it does not close all package dependencies or later-stage tools.

Archive identity is anchored to retained local metadata, not newly validated
signed repository metadata. No installation prefix transformations were
applied; the selected raw members match without them. Package archive bytes
are not source/license compatibility evidence or a security recommendation.
The separate base-runtime installation/restoration test, complete component
and rights review, terminal admission of job 22337, dedicated timing and
remaining scientific requirements are still unfinished.
