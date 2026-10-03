# Portable Presentation Rebuild

The accepted review can now be rebuilt in a relocated checkout using only
17 Git-managed files: the new script, accepted assembly inventory, original
main PDF and fourteen original figure PDFs. Add the existing 20,596-byte
fixed-length simulation figure without regenerating it. Its SHA256 is
`d4fc7ca1579234afff293672325b9b8e91c0fb9db43c41fdad3f1ed69c6ba4e2`.
The original assembly and rejected prototypes remain unchanged.

## Run

From the repository root, with PyMuPDF available:

```bash
python -B benchmark_tools/results/rebuild_publication_review.py --output /tmp/orthohmm-review-new
```

The directory must not exist. The output contains `document.pdf`,
`assembly.json` and appended-page preview PNGs. Input checksum checks finish
before output creation; existing directories are refused without writes.
Neither Git nor the original absolute source locations are required by the
rebuild itself. The accepted hash-pinned inventory supplies captions and PDF
hashes; historical producer manifests, raw results, large parameter metadata,
previous prototypes and their receipts are not loaded.

## Verification

```bash
python -B -m pytest -q benchmark_tools/results/test_rebuild_publication_review.py -p no:cacheprovider
```

Executed with Python 3.12.3, PyMuPDF 1.27.2.3 and pytest 9.1.1 in the existing
test environment: **8 passed in 10.20 seconds**. The tests export exactly the
17 staged Git objects into temporary relocated directories. Both normal and
optimized Python rebuilds match all 27 accepted pages in dimensions, text-word
geometry, pixels, bookmarks and link semantics. All six figure redirects and
45 restored file-URI actions remain present. Integrity/rendering checks are
explicit exceptions rather than optimization-removable assertions.

Other tests corrupt a PDF or the accepted inventory and check refusal before
output creation, and preserve an existing output marker. Each negative case
also runs under optimized Python. The earlier seven-case normal-mode panel
passed before adding optimized positive coverage; no failed test is hidden.

These are presentation tests using the same renderer, not independent
cross-viewer certification. Existing manual acceptance is reused through exact
page/link identity; no new manual review is asserted by the rebuild report.
External audit/data/code links retain their historical locations and still
need separate files. This does not reproduce plotting, inference, scoring or
statistical provenance, and is not a portable scientific release.

No timing/helper deployment, frozen settings, scientific scores, scheduling,
DGX, service or unrelated-job state changes. Controlled timing, final
scientific reconciliation and the complete release remain unfinished.
