# Documentation Dependency Security

The authenticated, read-only GitHub alert snapshot contains 21 open alerts:
one critical, seven high, eleven medium and two low. All target `docs/uv.lock`.
No claim that the inference runtime or installed host environments are free of
vulnerabilities follows from that repository-level inventory.

[Retained alert metadata](dependency_alerts_before_20260917.json) excludes
credentials, exploit descriptions and unrelated API fields. Its SHA256 is
`de817d17ecd97faae00abde6fb70f12fdfd2233143d21107b917b4f9726b150b`.
The installed `gh` is not the GitHub CLI; collection used the GitHub REST API
with the existing credential helper. Initial failed requests are not treated
as evidence of an empty alert inventory. The collector makes read-only calls
and does not dismiss alerts.

## Remediation Scope

| Package | Previous lock | Updated lock | Minimum constrained version |
| --- | --- | --- | --- |
| h11 | 0.14.0 | 0.16.0 | 0.16.0 |
| idna | 3.10 | 3.20 | 3.15 |
| Jinja2 | 3.1.4 | 3.1.6 | 3.1.6 |
| Pygments | 2.18.0 | 2.21.0 | 2.20.0 |
| requests | 2.32.3 | 2.34.2 | 2.33.0 |
| Starlette | 0.41.3 | 1.6.0 | 1.3.1 |
| urllib3 | 2.2.3 | 2.8.0 | 2.7.0 |

The first targeted resolution retained vulnerable Python3.9 branches. The
docs-only Python floor is now3.10 and `.python-version` selects3.12. OrthoHMM
package requirements, frozen native inference checkouts, benchmark environments
and the active DGX timing environment were not changed. CI pins uv0.12.15
and uses `--locked` rather than silently updating the dependency solution.
This follows uv's documented [universal resolution](https://docs.astral.sh/uv/concepts/resolution/)
and [locking behavior](https://docs.astral.sh/uv/concepts/projects/sync/).
Advisory-specific source links and affected ranges remain in the snapshot;
for example, the [Starlette advisory inventory](https://github.com/Kludex/starlette/security/advisories)
documents the HTTP-serving issues. A docs dependency is not automatically
reachable from OrthoHMM's inference CLI.

## Verification

- `audit_dependency_lock.py` checks every locked version, including alternate
  interpreter branches, against each of the21 retained advisory ranges.
  [Result](docs_dependency_range_audit_20260917.json): zero affected versions.
  This is not a complete vulnerability scanner or a reachability proof.
- Four focused tests cover alert-field allowlisting, package-name normalization,
  vulnerable alternate branches, patched versions and absent dependencies.
- A fresh isolated Python3.12.3 environment installed the locked docs/dev set.
  `uv pip check` found all34 installed packages compatible.
- Sphinx7.4.7 built all seven source documents, exit0. It reported14 diagnostics,
  including reStructuredText ERROR messages, duplicate targets, missing images
  and a language-setting warning in unchanged sources. This is not a clean
  documentation build. Raw warnings remain at
  `benchmarks/work/docs_security_build_warnings_v1.log`.
- The installed sphinx-autobuild served the built index over localhost with
  HTTP200. Its own subprocess was terminated after the smoke test. Raw log:
  `benchmarks/work/docs_preview_security_v1.log`. No remote preview was exposed.
- The updated lock SHA256 is
  `5b57557fc3b4a5dc1c8cf2d28cda0db2d67214d7da0abffb9cea7d7442cd5b31`.

The immediate [post-push snapshot](dependency_alerts_after_20260917.json) after
commitd29cdd8 still reports21open alerts. Local range checks do not establish
server-side closure; recheck after GitHub dependency processing. No alerts
were manually dismissed. Future advisories, existing documentation
errors, action-version hardening and historical environment exposure remain
separate release considerations. Do not silently update frozen environments
or call them safe because the current docs lock has been remediated.
