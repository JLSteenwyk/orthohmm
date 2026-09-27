# Repeatable Base Archive Staging

`benchmark_tools/stage_base_archives.py` replaces the manual archive-copy and
explicit-file construction step of the [base reconstruction](RECONSTRUCTED_BASE_FIXTURE_20260927.md).
It reads package identities from the committed reconstruction receipt, not
the original machine's archive or installed-metadata paths. Supply the 19
retained archives in a cache directory, using the filenames, HTTPS provider
URLs and hashes in that receipt. This tool does not download them.

```bash
python -m benchmark_tools.stage_base_archives \
  --receipt benchmark_tools/results/reconstructed_base_fixture_20260927.json \
  --cache /path/to/acquired-archives \
  --output /path/to/new-staging-directory
```

The destination must not exist. The complete cache is checked before output
creation; symlinked/missing archives, duplicate packages, malformed identities
and digest/size differences fail. Copied bytes are checked again. A changed
receipt fails before the explicit file or success receipt is written. Failed
directories remain available for inspection and are not automatically retried.
The tool uses only the Python standard library.

`explicit.txt` uses destination-specific `file:` URLs and MD5 fragments for
the documented offline Conda installation. SHA256, MD5 and byte counts are
all independently checked during staging. `staging.json` records relative
archive paths, package identities, the input-receipt SHA256 and explicit-file
identity. Restage from a relocated cache to regenerate valid file URLs;
simply moving `explicit.txt` does not make its URLs relocatable.

## Validation

All 22 focused tests pass, covering staging/relocation, paths with spaces,
invalid metadata, duplicate/empty package sets, corrupted/missing/symlinked
archives, changed receipts, output refusal and copy corruption. Two real
staging executions succeeded using all 19 archives (51,094,707 bytes):

| Local destination | Staging receipt SHA256 |
| --- | --- |
| `/tmp/orthohmm-base-staging-cli-20260927` | `2bce393aae020703611d61eef63f67607cec777c6ba0b35176492373d4eaffb4` |
| `/tmp/orthohmm-base-staging-relocated-20260927` | `6e730dd867b1eab0ad0e693d52fe310665634eb0b44d1772d00548f68bab8fe9` |

The second execution used the first output's archive directory as its cache.
All package/archive records agree between executions; explicit-file hashes
differ because their destination URLs differ. No original metadata path was
needed. These tests do not repeat Conda installation, the pip overlay or
scientific inference, and do not establish cross-host runtime compatibility.
The supplied receipt remains a trusted input, not an independently authenticated
provider signature. No security, dependency-closure or redistribution clearance
is implied. Job 22337 and all its pinned inputs were left unchanged.
