# Verified llvmlite Publishing Attestation

PyPI's [Integrity API](https://docs.pypi.org/api/integrity/) returned a
provenance bundle for the exact retained CPython 3.10 Linux x86_64 wheel.
[Raw query receipt](llvmlite_pypi_provenance_query_20260929.json) preserves
the response and local wheel identity. The statement's sole subject matches
both the retained filename and SHA256
`854941c2267fd4fc5b2ce02b8af8ecdffa79fb7784591d3a89370322039ea09f`.

Using the [documented verification tool](https://docs.pypi.org/attestations/consuming-attestations/),
`pypi-attestations` 0.0.30 with Sigstore 4.5.0 verified the local wheel against
that saved provenance and the expected repository `https://github.com/numba/llvmlite`.
The command exited 0 and reported OK. The
[verification receipt](llvmlite_attestation_verification_20260929.json) records
the exact command, output, wheel checksums and complete verifier package list.
The verifier was installed into a new temporary venv, not any inference,
reader or shared scientific environment. Verification used normal trust
refresh, not a disabled signature check or an offline bypass.

The [decoded identity](llvmlite_attestation_identity_20260929.json) binds
the subject check and selected certificate extensions to the verified receipt:

- Publishing run: `31511294686`, attempt 1.
- Workflow commit: `5d881fd159581b97cd8221029ed3b1e50203e32e`.
- Publisher: `numba/llvmlite`, `upload_packages.yml`, environment `pypi`.

The [public run metadata and immutable workflow](llvmlite_upload_run_20260929.json)
confirm a successful Upload Packages run. Its workflow gathers artifacts
from other runs and publishes them. It is not the LLVM compilation workflow.
This explains why the verified publishing commit differs from the earlier
wrapper release-tag commit; neither should silently replace the other.

## Boundary

This is a verified PyPI publish attestation with a null predicate, not SLSA
build provenance identifying all source and toolchain inputs. It closes the
exact wheel-to-publisher gap but does not yet bind the wheel to a specific
LLVM development package or build run. Next trace the upload's selected
artifact-producing runs and compare exact wheel digests where artifacts or
logs remain available. Missing records must remain explicit. No compiler,
scientific runtime, parameter, benchmark result or timing eligibility changed.
Publication and redistribution clearance remain unproven.
