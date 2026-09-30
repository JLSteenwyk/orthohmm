# Three Additional Historical QfO Context Layers

Follow-up to the [three release-tag inspection](TREEFAM_QFO_CONTAINER_INSPECTION_20260930.md).
This examines new tag targets, not a repeat download of those inspected layers.
The hypothesis remains that an untracked original tree/mapping input could
have entered a developer's Docker build context without appearing in Git.

## Actual Scope And Result

Used the unchanged inspection helper through its `inspect` function, selecting
the 2020 reference-generation update tag and two historical developer/test tags
from the public [Darwin registry listing](https://hub.docker.com/v2/repositories/qfobenchmark/darwin/tags?page_size=100).
Tag names do not authenticate data-release identity or build provenance.

| Tag | Copied-context compressed bytes | Tar members | NHX/mapping filename candidates |
| --- | ---: | ---: | ---: |
| `2020.2.1` | 79,538 | 69 | 0 |
| `silvia` | 96,099 | 80 | 0 |
| `test-sina` | 80,797 | 65 | 0 |

All three were accessible and had one recognized `COPY . /benchmark` layer.
Downloaded **256,434 bytes** total; all three compressed layers match their
manifest SHA256/size. Complete inventories contain **214 member entries**,
with no TreeFam-named member, NHX file or `treefam2reference.txt` candidate.
The layers contain benchmarking source/configuration. No container was run,
no archive member extracted, no person contacted and no family labels inferred.

Across the two inspections, six different copied-context layers now cover six
of the 40 listed tags: 501,755 compressed bytes and 425 tar member entries.
The other 34 tags, other layers/namespaces/registries, deleted images and
differently named or encoded contents were not exhaustively inspected.
This is not proof of global source unavailability. Original full TreeFam-A
trees/mapping and family-level uncertainty remain unresolved.

## Evidence And Reproduction

The [portable receipt](treefam_developer_container_context_receipt_20260930.json)
retains raw public tag/config/manifest JSON, descriptors, full inventories and
source pin. SHA256:
`eb35e4ae29056f7cc364c1769023afbac89349df2ef471377678fa7090e32201`.
Independent readback checks all 10 acquired evidence files, three layer
descriptors, and every complete inventory against freshly reread tar headers.
All 214 header records match exactly; no extraction was used.

Inspection source remains SHA256
`10e733cee7db9c21a28226cb9cee5b387c888e3e1d84eb292ee5bcd768d239a8`,
the tested/pushed version from `7adf23b9`. Its unchanged four component tests
were not repeated merely because a new tag list was selected.

```sh
/home/bizon/anaconda3/bin/python -B -c 'from pathlib import Path; from benchmark_tools.inspect_qfo_darwin_images import inspect; inspect(Path("benchmarks/work/treefam_developer_container_search_20260930"), ["2020.2.1", "silvia", "test-sina"])'
```

The destination must be fresh. Raw layer blobs remain local and uncommitted
under that directory. Mutable tag aliases can change; use the retained digests
to identify the inspected contents. The script's original bounded-response/
download rules remain unchanged. No scientific score, native pipeline,
prior receipt or running FAS source changed. Avoid repeating these layers
without a concrete new lead.
