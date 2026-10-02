# Private Historical Bootstrap and Frozen Base

The copied publication workflow can now acquire a fixed official Miniforge
installer, install it privately and use that bootstrap to reconstruct the
unchanged historical base. This removes the requirement to use the shared
Anaconda installation for this step. It is not a new scientific result,
controlled timing measurement, complete runtime closure or publication release.

## Source and Tests

Checkpoint `5016a6511e69fee092f235f4df47c247f4e88a5d` adds the stdlib bootstrap
controller and mocked tests. Corrected checkpoint
`36d31f80b2f085489afc6ca9d043f0435dbae63d` supplies the private bootstrap's
`bin` directory on PATH and extends validation/response guards. The offline
base controller similarly puts the supplied Conda entrypoint's directory first.
No scientific source, settings, scoring endpoint or dependency version changes.

The final four-module panel passes **228 tests in 20.33 seconds**, with zero
failures, errors or skips. It includes an actual tiny synthetic `env python`
launcher check; mocked cases do not execute an installer or Conda packages.
Earlier 217- and 227-case panels overlap this final panel. The retained final
JUnit file is 33,473 bytes, SHA-256
`97d533fc7fbf1c3406de48441ae8a888a86b32d4f3c2452d1e687fb18c4e002e`.
The expanded eight-module panel additionally passes **380 tests in 24.42
seconds**, zero failures/errors/skips, including both public acquisition
helpers, archive staging, source/handoff and scaling-controller tests. These
panels overlap; do not sum their counts. Times above are pytest CLI wall times;
the receipt separately records JUnit suite durations. The [final validation receipt](publication_private_bootstrap_validation_20261002.json)
binds both JUnit files and the execution/status receipts.

The corrected `native-wheels` export has 1,835 payloads: 43 frozen scientific
files and 1,792 workflow/support files, totaling 10,898,819 bytes. Its
742,317-byte external source index has SHA-256
`0c756ca4fb5dd3c1b31a357637f6ce1e24bb086496403d64dd315bc9cdcf016f`.
The copied tree verifies before and after execution outside the checkout with
isolated stdlib Python and Git absent from PATH; all 1,814 Python files compile
without executing them. The scientific revision remains
`7f3a9e40dd7e79f842cc2c11fb8b548f9a802806`.

## Actual Acquisition and Installation

The fixed provider is the [official Miniforge 25.3.1-0 release](https://github.com/conda-forge/miniforge/releases/tag/25.3.1-0),
not a latest-version endpoint. Both artifacts are actually downloaded once:

| Artifact | Bytes | SHA-256 |
| --- | ---: | --- |
| `Miniforge3-25.3.1-0-Linux-x86_64.sh` | 93,870,801 | `376b160ed8130820db0ab0f3826ac1fc85923647f75c1b8231166e3d559ab768` |
| Provider `.sha256` sidecar | 104 | `57be9d8415cd75326aff2a518bd6eda8d3a87ef28747256d51a8b77a47cc0a14` |

The exact sidecar agrees with the installer digest. HTTPS redirects are checked
against the fixed GitHub/release-asset hosts; temporary signed query strings
are not retained. The acquisition receipt is 1,801 bytes, SHA-256
`0ba81c03cf9aa32ab7267dd5ce27abe9c8f38c0726f0b272a5c01b1800e624da`.
This historical installer is not a current-security recommendation, and a
signed release commit does not authenticate every binary payload.

Installation uses a fresh private prefix, private HOME/cache/config, bounded
stages and [documented batch mode without default shell initialization](https://github.com/conda-forge/miniforge#unix-like-platforms-macos--linux).
Both installer and Conda-version stages complete successfully in the corrected
attempt. Conda reports **25.3.1**; 88 package metadata records have unique names.
Its 509-byte Conda entrypoint has SHA-256
`e0567dad4bcf85ac129422956f1faa41fc31b5648620d26805606397e337457b`.
Use the observed digest for this supplied entrypoint; do not assume another
installation's identity or treat the launcher as the full bootstrap closure.
The bootstrap completion receipt is 42,409 bytes, SHA-256
`91c2752539d263d9c7b30c289543b178027081e86b001bb0d4937a17d8bdfe8d`.

The newly private Conda then drives the copied offline base controller using
the already [publicly acquired 19 archives and pip wheel](PUBLICATION_BASE_ACQUISITION_20261002.md).
No base artifact is downloaded again. All four stages succeed once:

| Check | Result |
| --- | --- |
| Conda inventory | All 19 exact historical name/version/build triples |
| Interpreter | CPython 3.10.13, Linux x86-64, new private prefix |
| Pip overlay | Exactly pip 26.2.1; `pip check` passes |
| Pip payload | 475 wheel files match; generated RECORD excluded |
| Controller stages | Offline install, pip bootstrap, dependency check, snapshot; all zero |

The base completion receipt is 108,343 bytes, SHA-256
`0b6e43790ecd3ef5b6cd4caa8476a867633d3853e4c79cf1bb7d8687595606f0`.
Its 17,157,264-byte Python executable has SHA-256
`3f5a56d80012878d44e08cb5e1c9d314a405315933d6e87e96b197b4bfb34896`.
Installed executable hashes need not equal those of another prefix; metadata
equality is not whole-runtime byte identity or new scientific/native admission.

All five corrected export/verification/installation stages return zero. The
[machine receipt](publication_private_bootstrap_20261002.json), 75,610 bytes,
SHA-256 `4e881c067b5c19630fed8592d04cc8a62cb65bf21bddebe3fa9b5fe821900970`,
binds their commands/logs, downloads, inventories, inputs, failures and source
anchor. Reader commands are in the [source-component guide](../PUBLICATION_SOURCE_COMPONENT.md#private-historical-conda-bootstrap).

## Failures Retained

The first receipt checker expects a bare filename but the exact sidecar uses
`./filename`; it stops after successful acquisition, before installation.
Correct only readback and reuse both verified downloads.

The first real installation finishes, but Conda's `#!/usr/bin/env python`
launcher cannot find Python on the original restricted PATH. Version checking
fails, so no bootstrap success is emitted. Preserve the failed prefix and logs.
After the tested PATH correction, make one justified fresh-prefix installation;
do not overwrite the failed attempt or silently count it as successful.

Two subsequent one-off receipt checks fail: one includes pip release metadata
among base archives, and one expects full base fields in the CLI's summary
instead of its digest-bound `complete.json`. Correct only role selection and
readback. Neither successful bootstrap nor successful base installation is
restarted. The final receipt binds all five failure records covering these
four distinct issues, including the controller's separate PATH-failure record.

## CI and Remaining Scope

At 17:20:38 UTC, the preceding pushed revision
`27c0b138c278a8c2abef5ad15ec856c5dd07154a` has completed
[CI run 37037578636](https://github.com/JLSteenwyk/orthohmm/actions/runs/37037578636)
successfully: all eight observed jobs are successful, including four fast
Python-version jobs, coverage, docs, wheel and Linux native diagnostics.
The [status receipt](ci_private_bootstrap_prior_revision_20261002.json) records
exact IDs/statuses; no remote test counts were extracted. That run validates
the preceding CI portability fix, not these newer bootstrap changes.

No shared environment, service, unrelated process, DGX, quiet-window question
or contention poll is involved. Completed full OrthoBench/native fixtures and
collector admissions are reused, not restarted. Production timing remains
deferred; the goal remains active.

Installer identity and metadata checks are not complete installed-bootstrap,
prefix-transformed base, native/OS, source attribution, security or rights
closure. Same-host copying is not cross-host restoration or a syscall sandbox.
Unpublished project-wheel/tool distribution, remaining QfO uncertainty,
comparable resources, final manuscript/archive/release/DOI and overall
publication readiness remain unfinished.
