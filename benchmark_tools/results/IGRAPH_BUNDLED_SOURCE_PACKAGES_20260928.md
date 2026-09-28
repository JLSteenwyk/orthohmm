# Signed Distribution Source Packages Recovered

The version-specific AlmaLinux 8.10 vault contains all three source RPMs
declared by the binary packages matched to the retained igraph wheel.
This resolves the preceding candidate-location failures, not by substituting
newer releases but by locating the named source packages.

| Source RPM | Download bytes | Regular members | Patch files |
| --- | ---: | ---: | ---: |
| libxml2-2.9.7-21.el8_10.3.src.rpm | 5,492,177 | 37 | 35 |
| xz-5.2.4-4.el8_6.src.rpm | 1,075,981 | 5 | 1 |
| gcc-8.5.0-28.el8_10.alma.1.src.rpm | 65,715,094 | 63 | 59 |

The [source-material inventory](igraph_bundled_source_material_20260928.json)
binds download URLs, RPM metadata, every member and eleven selected notice/header
exports. A [second extraction](igraph_bundled_source_recheck_20260928.json)
reproduced all 105 member payloads; source RPM names also match the binary
packages' `SOURCERPM` fields. Extraction rejects nonregular, duplicate or
non-flat source-package members. No specification, patch or package script was
executed and no compiler or runtime was installed or changed.

Raw packages, nested upstream archives, specifications, patches and selected
exports remain under `benchmarks/work/igraph_bundled_sources_20260928`.
The repository receives inventories and analysis, not large source archives.

## Signature Verification

All six source/binary RPMs report valid header/payload digests and RSA/SHA256
signatures under an isolated temporary RPM database. The
[signature receipt](igraph_bundled_rpm_signatures_20260928.json) retains the
complete output, key bytes/hash, primary fingerprints and tool identities.
No system RPM database or GPG keyring was modified.

The two primary fingerprints were checked against the
[official AlmaLinux security page](https://almalinux.org/security/#gpg-keys):
`5E9B8F5617B5066CE92057C3488FCF7C3ABB34F8` and
`BC5EDDCADF502C077F1582882AE81E8ACED7258B`. The older XZ packages use the earlier
signing subkey; the newer packages use the replacement key described in
[AlmaLinux's key-transition notice](https://almalinux.org/blog/2023-12-20-almalinux-8-key-update/).
This is an official-HTTPS trust basis, not an offline trust ceremony. The
security page was readable through web retrieval; a separate local raw-HTML
download returned 403 and was not bypassed or represented as archived evidence.

The initial check against the machine's unconfigured RPM trust database returned
nonzero; the isolated, explicitly pinned-key checks above are the retained
verification, not an assertion that an untrusted key error proves corruption.

## Notice Material And Remaining Limits

Source libxml2 `Copyright` and XZ `COPYING` are byte-identical to the two notice
texts recovered from the binary RPMs. The GCC archive includes `COPYING3` and
`COPYING.RUNTIME`; its `libgomp/libgomp.h` explicitly identifies GPL version 3
or later with GCC Runtime Library Exception 3.1 and refers to those two texts.
Their exact bytes are retained. This resolves the missing source-text lead for
libgomp, not an automatic compatibility determination for the whole wheel.

The distributions contain many patches; their presence must not be confused
with proof that each patch applies to this particular library or was compiled
into the wheel. Neither a filename match nor valid signatures reproduce a
build. The prior selected-section equality remains valid, but full
source-to-binary reproduction, complete transitive attribution and final release
compatibility review are still open. No security remediation or publication
readiness is claimed, and frozen scientific environments remain unchanged.
