# Patched Private Timing Candidate

Built a separate candidate at
`benchmarks/work/threadripper_patched_runtime_20260928/venv` using the reconstructed
private Python 3.10.13 base. Neither the shared installation nor the preceding
private candidate was changed. This implements the deployment revision proposed
in [the security review](PRIVATE_TIMING_SECURITY_REVIEW_20260928.md).

## Installation

The [installation receipt](threadripper_patched_runtime_20260928.json) records
the source, baseline checksum, downloaded wheels, commands, logs and explicit
deployment changes. The [full hash lock](threadripper_patched_requirements_20260928.txt)
contains all 30 distributions. Downloads use the explicit PyPI index and binary
wheels only; installation uses the local wheels with no index, no dependency
resolution and required hashes. All five stages succeeded: download, venv
creation, installation, pip check and the declared native import probe.

Relative to the preceding aligned candidate, only Pillow (12.2.0 to 12.3.0)
and setuptools (81.0.0 to 83.0.0) version selections change. Packaging stays
26.1, the version previously aligned with historical payload bytes. The receipt
also records its difference from the older baseline's stale 26.0 metadata.
All other selected scientific versions and the frozen scientific source remain
unchanged. Live metadata matches the full 30-package selection.

## Bounded Validation

The [fresh repository snapshot](dependency_alerts_patched_timing_20260928.json)
contains 40 open alerts. Its 14 additional alerts name the v2 historical lock;
they do not introduce a new advisory identity. No historical alerts were
dismissed or locks rewritten. The [live range check](patched_timing_advisory_ranges_20260928.json)
finds zero matching advisory ranges in the patched candidate. Absent packages
are distinguished from installed versions outside a reported range. This is
not a comprehensive security, reachability, native-library or OS audit.

The [import comparison](threadripper_patched_import_comparison_20260928.json)
rechecks all 40 frozen core files and compares 710 historical imported files:
704 have identical hashes. Changed imports are `PIL.Image`, `PIL.ImageFile`,
`PIL.PngImagePlugin`, `PIL._imaging`, `PIL._version`, and `yaml._yaml`.
Six unrelated editable startup finders remain omitted. The import probe sees
no shared Conda prefix in Python paths or imported module paths; this is not
an ELF dependency closure or file-access sandbox. Setuptools is inventory-
checked but this probe does not establish its complete payload behavior.

Thirty focused tests pass across environment building, advisory auditing and
import comparison. The new builder flag is explicit and leaves the historical
selection mode unchanged; tests verify that only the three documented package
selections differ from historical metadata and that baseline data is not mutated.

## Remaining Admission Work

This candidate is installed, not admitted for controlled timing. Validate native
fixture outputs and collector v5 against this deployment, document the YAML
build difference, then update runtime binding and the execution amendment.
Final resource accounting and whole-run quiet-host requirements still apply.
No production timing or scientific benchmark was rerun, no superiority claim
changed, and no DGX access occurred. The candidate depends on its private base
prefix and host OS; it is not a self-contained archive.
