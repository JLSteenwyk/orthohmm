# Bundled igraph Libraries: Distribution Evidence

Inspected the exact retained igraph 1.0.0 x86-64 wheel without importing it,
executing its libraries or modifying any benchmark environment. Its ELF debug
links identify more specific distribution builds than the shared-library names.
Downloaded matching binary RPM candidates from AlmaLinux's HTTPS repository.

| Wheel component | RPM candidate | Source RPM declared by candidate | Embedded notice recovered |
| --- | --- | --- | --- |
| libxml2 | libxml2-2.9.7-21.el8_10.3.x86_64 | libxml2-2.9.7-21.el8_10.3.src.rpm | Copyright, 1,289 bytes |
| liblzma | xz-libs-5.2.4-4.el8_6.x86_64 | xz-5.2.4-4.el8_6.src.rpm | COPYING, 2,775 bytes |
| libgomp | libgomp-8.5.0-28.el8_10.alma.1.x86_64 | gcc-8.5.0-28.el8_10.alma.1.src.rpm | No standalone license member in the listed binary package |

For all three libraries, `.text`, `.rodata`, `.eh_frame`, `.gnu_debuglink`
and `.note.gnu.build-id` match the respective RPM library byte-for-byte.
All three whole-object hashes differ. This is stronger component attribution
than a version-looking filename, but not complete binary equality or a source
build attestation. No conclusion is drawn about which unexamined sections differ.

The [section comparison](igraph_bundled_rpm_sections_20260928.json) records
wheel and RPM identities, direct download URLs, member hashes, all fifteen
section comparisons, tool hashes, RPM metadata and notice-file identities.
An [independent second extraction](igraph_bundled_rpm_recheck_20260928.json)
reproduced all thirty section payload hashes and both notice texts.
Raw RPMs and notice copies stay under
`benchmarks/work/igraph_bundled_rpms_20260928`, not in the source repository.

## Public Build Context

The annotated upstream 1.0.0 tag resolves to commit
`b16f27618674dd1913007a52855b76075802cbf9`. Its
[build workflow](https://github.com/igraph/python-igraph/blob/b16f27618674dd1913007a52855b76075802cbf9/.github/workflows/build.yml)
uses cibuildwheel 3.2.1 and installs distribution development packages, including
libxml2, before building the C core. This is consistent with a distribution
dependency origin; it does not establish which job produced the retained wheel
or pin the container digest. The downloaded tag response and workflow are
hash-bound in the [context receipt](igraph_bundled_rpm_context_20260928.json).

RPM `LICENSE` metadata declares MIT for libxml2, Public Domain for xz-libs,
and a multi-license expression with exceptions for libgomp. These declarations
are retained as provider metadata, not legal compatibility determinations.
No RPM was installed and no package script was run. Package signatures were not
verified against an independently pinned trust root, so HTTPS acquisition and
checksum reproduction are not a distribution-signature attestation.

Six candidate source-RPM URLs under AlmaLinux BaseOS/AppStream returned 404.
Those failures are recorded; they do not establish global source unavailability.
Corresponding source acquisition, patches/build-chain verification, missing
libgomp notice attribution and broader compiled/transitive component review
remain open. The wheel is not cleared for redistribution by this analysis.

## Timing Remains Separate

A [fresh local process capture](threadripper_timing_recheck_after_stream_20260928.json)
measured 92.5922 competing CPU-core equivalents despite an empty Slurm queue.
It retained thirteen sampling errors and three unmatched foreign processes.
The largest observed service groups were the BAli-Phy independent chains,
matched-interval analyses and Neocallimastix guide analyses, each near sixteen
cores. This is incomplete sampled demand, not a complete census or permission
to interrupt them. No production timing or DGX work was started.
