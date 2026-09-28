# Reader Dependency Security Upgrade

The [September 28 environment-scope review](RELEASE_SECURITY_SCOPE_20260928.md)
rechecks the live reader metadata and distinguishes this patched installation
from the benchmark interpreters, which still match reported advisory ranges.
The reader result is not security clearance for those separate environments.

The repository's September 27 advisory snapshot increased from 11 to 12 open
alerts after the reader lock was added. The new alert is CVE-2025-68463 /
GHSA-x3vf-39hj-gxr4, concerning untrusted XML processed by Bio.Entrez. The
snapshot reports affected versions `<= 1.86` but no first-patched version.
The upstream [Biopython 1.87 release notes](https://raw.githubusercontent.com/biopython/biopython/biopython-187/NEWS.rst)
explicitly identify the fix. The [upstream issue](https://github.com/biopython/biopython/issues/5109)
describes XML-driven external requests. No exploit or reachability test was run.

## Executed Validation

A new [reader-only lock](publication_reader_requirements_20260927_v2.txt)
replaces Biopython 1.86 with 1.87. Other four versions are unchanged. The
CPython 3.10 Linux wheel was downloaded from PyPI, then installed offline,
binary-only and hash-required in `/tmp/orthohmm-reader-patched-20260927/venv`.
Neither the historical reader installation nor the frozen inference changed.

The [validation receipt](reader_upgrade_validation_20260927.json) records
five installed distributions, wheel and lock identities, pip check, before/after
payload audits, all five report comparisons, execution command and file trace.
All four independent readers pass on the retained relocated 16-gene native
fixture. Scientific report fields are identical. The recorded Biopython
version changes from 1.86 to 1.87; internal report references are rehashed
before being matched by report identity. External checked records are rehashed
and must remain identical. All 2,702 audited package payload files match their
wheels both before and after readback. Pip check reports no broken requirements.
Thirty-five focused tests pass across the upgrade comparator, advisory
inventory, wheel-payload auditor and source exporter, including 14 new tests.
Eighteen top-level receipt file identities were rechecked after validation.

The first comparator attempt failed on the expected Biopython version and
relocated internal report references, after the readers themselves succeeded.
Its outputs remain under the fresh runtime's initial `readback/` directory.
The corrected comparator has explicit version reporting and validated pointer
normalization, with scientific equality unchanged. Successful execution is
under `validation_v2/`; no original output was overwritten. This was a reader
validation repeat, not a native inference or benchmark-scoring retry.

All 12 retained repository advisory ranges were evaluated against the verified
live installation: zero affected versions. Historical locks and their GitHub
alerts remain intact. This is not an assertion that all alerts are closed,
that all dependencies have been screened, or that historical environments
are safe for untrusted input.

## Reproduction

After installing the new lock in a fresh environment with its five verified
wheels, run the repository's validator with absolute paths. It refuses an
existing output directory or receipt:

```bash
python -m benchmark_tools.validate_reader_upgrade \
  --runtime /fresh/reader \
  --readers /exported/independent-readers \
  --native /retained/native \
  --previous /previous/readback \
  --snapshot /retained/repository-alerts.json \
  --lock benchmark_tools/results/publication_reader_requirements_20260927_v2.txt \
  --output-directory /fresh/validation \
  --output /fresh/validation-receipt.json
```

The same-host fixture has no satellite merges. This result does not establish
full-dataset version equivalence, independent accuracy, OS/native-library
security, cross-host portability, or publication readiness. Payload auditing
excludes generated RECORD/bytecode and relocated non-site wheel data. The
trace is evidence, not an access sandbox. Combined release packaging and the
previously listed scientific/timing gates remain open.
