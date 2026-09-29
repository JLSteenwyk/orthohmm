# Candidate Release PyPI Advisory Check

The [release-specific snapshot](release_pypi_advisories_20260928.json) expands
the earlier check of repository-reported advisory ranges to every package in
the candidate recovery and patched-reader locks. It uses the
[PyPI release JSON API](https://docs.pypi.org/api/json/#get-a-release), not the
latest-project endpoint. No installation, upgrade or scientific run occurred.

| Lock | Pinned packages | Public releases queried | Unresolved | Active advisories returned |
|---|---:|---:|---:|---:|
| Recovery environment | 11 | 10 | 1 | 0 for the 10 queried releases |
| Patched reader | 5 | 5 | 0 | 0 |

Across both locks there are 12 distinct name/version pairs, of which 11 have
successful public responses. All locked hashes for those public releases match
published artifacts. No known advisory was returned for them. This is a
time-specific advisory observation, not evidence that vulnerabilities are absent.

The `orthohmm==0.5.0` endpoint returned HTTP 404. The local source-built wheel
therefore remains unresolved by this public-package check; it is not counted
as a clean release or equated with any differently versioned public package.
Existing source and wheel provenance remain separate evidence.

The [compressed response bundle](release_pypi_snapshot_20260928.json.gz)
retains the exact successful JSON response bodies and original report, indexed
by filename. The result records retrieval times, URLs, response hashes, lock
hashes and auditor identity. All 18 referenced file records were rechecked;
the compressed bundle was round-tripped against the source files. Network
errors are explicitly unresolved, withdrawn advisories remain in the record,
and malformed/mismatched responses cannot become empty clean results.

Thirty-eight focused tests pass across this auditor, the reader-upgrade
validator and the earlier advisory comparator. An initial test command named
a nonexistent test file and collected no tests; the corrected command passed.
The initial API snapshot is retained locally; the final v2 snapshot binds the
auditor after removal of an unused import and expanded tests, with unchanged
substantive findings.

This does not resolve transitive dependencies, inspect live package payloads,
assess exploitability, scan embedded native libraries or the OS, establish
bootstrap trust, or validate full-dataset equivalence. Historical benchmark
environments and GitHub alerts were not modified. Release readiness remains
unproven.

```bash
python -m benchmark_tools.audit_pypi_releases \
  --lock benchmark_tools/results/publication_recovery_requirements_20260926.txt \
  --lock benchmark_tools/results/publication_reader_requirements_20260927_v2.txt \
  --output-directory /fresh/path/pypi-audit
```
