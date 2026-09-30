# Selected QfO Container Build-Context Inspection

## Question And Scope

Could original TreeFam-A release-7 trees or `treefam2reference.txt` have entered
a public QfO image through an untracked local build context, despite their
absence from the previously inspected Git history?

The retained [QfO repository](https://github.com/qfo/benchmark-webservice) at
`c0854a96c1a0fd7f2a891d971af0863002fabc90` configures the Darwin image tag
`2022.1`. Its Dockerfile copies the build context into `/benchmark`;
`.dockerignore` excludes `reference_data` but does not explicitly exclude
`data/treefam`. This motivated inspecting actual image layers rather than
inferring absence from Git or the Docker recipe.

Used anonymous, read-only public registry access and the
[Registry API](https://distribution.github.io/distribution/spec/api/).
No credentials were requested or retained. No Docker daemon, container runtime,
image execution, archive-member extraction or service change was used.

## Actual Result

The public [Darwin tag listing](https://hub.docker.com/v2/repositories/qfobenchmark/darwin/tags?page_size=100)
returned 40 tags in one page. Selected **three**: early `2020.1`, matching
reference-generation `2020.2`, and configuration-selected `2022.1`.
Tag names alone do not establish which reference files were used historically.

| Tag | Copied-context compressed bytes | Tar members | NHX/mapping filename candidates |
| --- | ---: | ---: | ---: |
| `2020.1` | 71,609 | 64 | 0 |
| `2020.2` | 87,520 | 72 | 0 |
| `2022.1` | 86,192 | 75 | 0 |

All three selected manifest/config histories contain one recognizable
`COPY ... /benchmark` layer. Empty history entries were excluded when mapping
history to manifest layers. The OCI index for `2022.1` was resolved to its
unique Linux/amd64 image, with descriptor size and SHA256 verified.
The other two tags directly provide image manifests.

Downloaded only the three copied-context layers: **245,321 bytes** total.
Each compressed blob matches its manifest size and SHA256. Complete inventories
retain **211 member entries**, including directories and whiteouts. These layers
contain benchmarking source/configuration, not TreeFam-named files, NHX files,
or `treefam2reference.txt`. The inventory includes the reference-generation
script `generateData/AddReconciledTree.drw`, not its missing original inputs.

This is a negative finding for those three **copied-context layers only**.
The remaining 37 listed tags, other image layers, Python/FAS images, other
registries, differently named/encoded content and deleted historical images
were not exhaustively searched. It is not proof that no public source copy exists.
No original family labels or uncertainty estimates are enabled by this result.

## Reproduction And Evidence

Executed [inspection helper](../inspect_qfo_darwin_images.py):

```bash
/home/bizon/anaconda3/bin/python -m benchmark_tools.inspect_qfo_darwin_images \
  benchmarks/work/treefam_container_search_20260930
```

The destination must be new; existing attempts are not overwritten.
Prospective download limits: 2 MB per metadata response, 32 MB per selected
context layer, and 96 MB total unique context blobs. Oversized or unrecognized
artifacts remain unresolved rather than becoming absence findings. Mutable tags
can change; reproduce the recorded contents using the retained manifest/blob
digests rather than assuming today's tag aliases are immutable.

Executed source SHA256:
`10e733cee7db9c21a28226cb9cee5b387c888e3e1d84eb292ee5bcd768d239a8`.

The [portable receipt](treefam_container_context_receipt_20260930.json) includes
the raw public tag/config/manifest JSON, full inventories, selected descriptors,
source pin and independent file-size/checksum readback. Its SHA256 is
`762a6d7cfa9c096f165b1750b673dfa8b5f299fb2fc0348f46d9c2bdefcb3f42`.
All 11 downloaded registry evidence files and all three inventories were read
back; selected layer hashes/sizes and all inventory counts match. An initial
readback command used `json` rather than `.json` as its suffix filter and failed
before writing a receipt; correcting that command produced the retained result
without repeating network acquisition or changing the inspection helper.

Four focused tests pass: digest/size rejection, exact nonempty history mapping,
history mismatch rejection and inventory-only handling of unsafe paths/links.
These tests do not constitute scientific scoring or runtime admission.

Raw blobs and metadata remain local under the isolated work directory above;
raw container layers are not committed. No scientific defaults, benchmark
scores, original source pins or prior search artifacts changed. Family-level
TreeFam uncertainty remains unavailable. Controlled timing remains deferred
under the user's latest direction; no contention probe or scheduling question
was repeated, and unrelated analyses were left alone.
