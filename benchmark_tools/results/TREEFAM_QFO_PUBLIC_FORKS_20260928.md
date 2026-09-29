# Public QfO Fork Source Search

Extended the earlier upstream-history inspection to the public forks returned
by GitHub's unauthenticated forks API. All eight returned repositories had
complete, nontruncated recursive default-branch inventories (126-194 entries).
The response contained fewer than the requested 100 repositories, so this
single listing page was exhausted. This is not a search of other branches,
full fork histories, deleted refs or private repositories.

Five forks contain `reference_data/data/TreeFam-A.json`, all at Git blob
`27849bc6a3e36d5403f70d53dfff5cd70ecf1257`. The other three inventories have
no filenames matching TreeFam or NHX. Downloaded the unique candidate blob,
decoded it and independently verified Git's length-prefixed SHA-1 identity.
Its 284 bytes define an aggregation visualization with TPR/PPV axes and an
empty participant list: no tree, family membership or mapping data.

The [search receipt](treefam_qfo_forks_20260928.json) records repositories,
default branches, tree IDs, candidate paths, API URLs, retrieval times and
SHA-256 identities. The [compressed responses](treefam_qfo_fork_responses_20260928.json.gz)
retain eleven exact response/content strings: the fork list, eight recursive
trees, candidate blob response and decoded metadata. Every source hash was
rechecked and the compressed bundle round-tripped. Four focused tests pass.

Sources: [public fork listing](https://api.github.com/repos/qfo/benchmark-webservice/forks?per_page=100),
[candidate blob](https://api.github.com/repos/inab/benchmark-webservice-refdataset-fetch-issue/git/blobs/27849bc6a3e36d5403f70d53dfff5cd70ecf1257),
and [upstream repository](https://github.com/qfo/benchmark-webservice).
The browser tool could not open the API listing; the normal HTTPS API client
succeeded. No access controls were bypassed or credentials used.

Neither original TreeFam-A7 trees nor `treefam2reference.txt` were recovered.
No family labels, uncertainty estimates or benchmark scores changed, and no
one was contacted. Filename-negative inventories do not establish absence of
data under other names or in other archives. The missing-source limitation
remains in force.

```sh
python -B -m benchmark_tools.inspect_qfo_public_forks --output /fresh/path/qfo-forks
```

This repeats the public inventory query, not historical-byte reconstruction.
Future fork/default-branch inventories may change. The candidate download and
Git-blob check are recorded separately in the retained receipt.
