# Public Acquisition of the Frozen Publication Base

Status: 2 October 2026. The exact historical base packages and pip bootstrap
wheel can now be acquired using a copied source component and public HTTPS
providers. Actual fresh downloads pass the existing offline installer's
preflight. No installation, native inference, benchmark or timing was rerun;
the historical runtime is not promoted as security-cleared or current.

## Executed Workflow

Committed source `254ab20bcf0306d109df725d83ec414c81917bd7` adds
[the acquisition controller](../acquire_publication_base.py), tests and
[reader commands](../PUBLICATION_SOURCE_COMPONENT.md#offline-historical-base-installation).
It consumes the exact frozen reconstruction receipt, not mutable local Conda
metadata, latest-version solving or historical archive paths. All 19 package
URLs, sizes and hashes come from that authenticated receipt. For pip it uses
the [documented release-specific PyPI JSON API](https://docs.pypi.org/api/json/#get-a-release)
and requires the one recorded wheel name, size and SHA-256, not an alternative.

The controller rejects non-HTTPS, credentials, query/fragment URLs, unexpected
hosts/ports and unsafe redirects. Every response is bounded, declared lengths
are checked when present, and all binary downloads must match retained hashes
and sizes before their partial files are admitted. Archive MD5 is an additional
identity check, not an authenticity argument. Failure preserves partial files
and records without completion, retry, installation or overwriting an old path.

The [executed machine receipt](publication_base_acquisition_20261002.json) is
47,225 bytes, SHA-256
`c6c912a7f085c2d12e7154d922079874a5e6e92a30b045b9bf4a638dbf0461ac`.
All five actual stages return zero: committed export, copied verification,
public acquisition, offline-input preflight and final copied verification.

| Check | Result |
| --- | --- |
| Copied source | 1,828 payloads: 43 frozen scientific and 1,785 workflow files |
| Syntax compilation | 1,810 Python files, without importing scientific source |
| Actual acquisition | 19 exact base archives and one exact pip wheel; 52,911,339 artifact bytes |
| Provider requests | 21: 19 archives, release metadata and wheel |
| Providers observed | `repo.anaconda.com`, `conda.anaconda.org`, `pypi.org`, `files.pythonhosted.org` |
| Installer preflight | All 19 packages, pip wheel and supplied trusted Conda entrypoint accepted; no output prefix created |
| Source after acquisition | Entire externally anchored copied component unchanged |

Execution is from a fresh copied tree outside the repository with isolated
standard-library Python. Component verification has no Git on PATH. Acquired
artifacts and logs remain in the receipt's private destination, outside the
immutable component; no downloaded archive/wheel is committed or publicly
redistributed. Source index: 739,167 bytes,
SHA-256 `8961e13262b4548cd24d7b5d9191fe6dda509edffcb28901cbe82dc838777d26`.
Acquisition's own `complete.json`: 25,906 bytes,
SHA-256 `fbdeb529cec2844967985aaf8d473bc93e6c65c5a28afcd88837b9ba871db9e4`.

Pip release metadata is retained separately, 6,683 bytes,
SHA-256 `4f5abeb5acacd5eaaad5e6982f5721dee3592795e59fc5a4067fcca92eed70e1`.
Its selected wheel is 1,816,632 bytes,
SHA-256 `71138adf1f4ca900cdb7d289c21b7494329f2332b6d85f0e1c42108c0384ed3e`.
Future metadata updates cannot authorize a different wheel or alter the
frozen package inventory.

## Validation and Remaining Scope

The final acquisition/installer/stager/source/handoff test panel passes
**214 tests in 14.22s**, zero failures/errors/skips. Its 29,876-byte XML is
pinned in the machine receipt at SHA-256
`17c909dc39158016c0334e55b274c4b3275849a1cd04f4bbb9f42f6144c564df`.
The preceding 92-case panel overlaps; fixture downloads/installers are
synthetic and distinct from the actual public-provider execution.

This closes the manual base-artifact acquisition step, not whole native-study
distribution. Conda bootstrap and transitive bootstrap dependencies, scientific
wheel/tool sets, separately acquired raw datasets, shared OS libraries and
complete source/rights obligations remain outside this downloader. Hashes and
HTTPS do not prove supply-chain signatures, security or redistribution rights.
This is same-host relocation, not a cross-host or syscall-sandbox test; its
timeout bounds socket operations rather than total workflow duration.

Reuse [completed installation and scientific fixture evidence](OFFLINE_BASE_CONTROLLER_20261002.md)
and the admitted full OrthoBench archive reproduction. No scientific settings,
scores, endpoints or uncertainty estimates changed. Controlled Threadripper
timing, remaining QfO uncertainty, complete runtime/rights closure, final
manuscript reconciliation and versioned public deposition/DOI remain open.
There was no host-contention poll, new quiet-window question, DGX operation or
unrelated process/service change. The full publication goal remains active.
