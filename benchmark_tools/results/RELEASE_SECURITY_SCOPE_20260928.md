# Release Security Scope

The later [private timing review](PRIVATE_TIMING_SECURITY_REVIEW_20260928.md)
supersedes the repository-wide alert count below: its snapshot has 26 open
alerts. The new private candidate matches 14 unique advisory ranges, including
13 Pillow advisories absent from the older snapshot. Earlier environment
comparisons below remain historical, bounded checks, not current clearance.

The subsequent [release-specific PyPI check](RELEASE_PYPI_ADVISORIES_20260928.md)
queries every pin in the candidate recovery and reader locks: 11 distinct
public releases return no known advisories and matching artifact hashes;
the local OrthoHMM 0.5.0 release remains unresolved by PyPI (HTTP 404).
This broadens package-advisory coverage, not native-library or OS clearance.

The [fresh repository snapshot](dependency_alerts_release_review_20260928.json)
contains 12 open alerts, all on two historical requirement files. No alert
was dismissed and neither historical lock was edited. This snapshot is a
bounded advisory inventory, not a comprehensive dependency vulnerability scan.

Live distribution metadata was collected without installing or upgrading
packages. The [environment comparison](release_advisory_environment_matrix_20260928.json)
records exact interpreter identities, commands, inventories and advisory-range
comparisons. The [recovery recheck](recovery_advisory_ranges_20260928.json)
also verifies its live inventory against the retained installation report.

| Environment | pip | setuptools | Biopython | Matching advisory ranges | Alerted package absent |
|---|---|---|---|---:|---|
| Benchmark OrthoHMM interpreter | 26.0.1 | 81.0.0 | 1.86 | 6 | None |
| Benchmark OrthoFinder interpreter | 24.0 | Not installed | 1.87 | 6 | setuptools |
| Patched independent reader | 26.2.1 | 83.0.0 | 1.87 | 0 | None |
| Recovery installation | 26.2.1 | 83.0.0 | Not installed | 0 | Biopython |

Counts are advisory-version matches, not confirmed exploitable paths. Multiple
advisories can affect one package. Absent packages are reported separately,
not represented as patched installations. Reader metadata matches its prior
validated five-package inventory; package-payload checks and scientific
fixtures were not rerun in this recheck.

The primary advisories have distinct scopes. The
[Biopython advisory](https://github.com/advisories/GHSA-x3vf-39hj-gxr4)
concerns XML handling in Bio.Entrez; the
[pip advisory](https://github.com/advisories/GHSA-qwm4-qh6w-59xr)
concerns package URL handling; the
[setuptools advisory](https://github.com/pypa/setuptools/security/advisories/GHSA-h35f-9h28-mq5c)
concerns source-distribution exclusions and Unicode normalization on macOS.
A version-range match alone does not establish applicability to a Linux
inference invocation. No exploit, reachability or platform-impact assessment
was performed here.

## Release Requirements

- Do not present the frozen benchmark interpreters or historical locks as
  patched, general-purpose release environments or use them for new installs.
- Keep separate the benchmark provenance, the
  [patched reader installation](READER_SECURITY_UPGRADE_20260927.md), and the
  validated recovery installation. This audit does not merge their evidence
  or establish full-dataset equivalence after dependency changes.
- Any remediation of a frozen scientific runtime needs a new pinned
  environment and appropriate equivalence checks, not an in-place upgrade.
- Retain the original alerts and limitations. New release packaging still
  needs comprehensive dependency, compiled-library, OS and bootstrap-chain
  review; zero matches in this small advisory set is not security clearance.

The existing recovery advisory tests pass (five tests). No global environment,
native scientific runtime, requirement lock or benchmark result was changed.
