# Observed Base Runtime Package Attribution

The [local package attribution](integrated_base_runtime_packages_20260927.json)
extends the [search-stage snapshot](INTEGRATED_SEARCH_LIBRARIES_20260927.md)
to its 46 non-wheel library paths and Python executable. Each of these 47
files has exactly one local package owner. This is a bounded inventory of
observed files, not the full runtime or a new installation recommendation.

| Local provider | Version | Build |
| --- | --- | --- |
| Conda python | 3.10.13 | h955ad1f_0 |
| Conda bzip2 | 1.0.8 | h7b6447c_0 |
| Conda libffi | 3.4.4 | h6a678d5_0 |
| Conda libgcc-ng | 11.2.0 | h1234567_1 |
| Conda libgomp | 11.2.0 | h1234567_1 |
| Conda libstdcxx-ng | 13.2.0 | hc0a3c3a_7 |
| Conda libuuid | 1.41.5 | h5eee18b_0 |
| Conda openssl | 3.0.16 | h5eee18b_0 |
| Conda xz | 5.6.4 | h5eee18b_1 |
| Conda zlib | 1.2.13 | h5eee18b_0 |
| Ubuntu libc6:amd64 | 2.39-0ubuntu8.4 | Not a Conda build |

## Checks Performed

Rehashed every observed file against the earlier snapshot. Conda ownership
uses exact relative paths in `/home/bizon/anaconda3/conda-meta/*.json` file
lists, not a basename match. The receipt retains each owner record's identity,
name/version/build, declared license, package archive digests and applicable
installed-file entry. When present, `sha256_in_prefix` was compared directly;
otherwise an unmodified hardlink entry without a prefix placeholder may use
its recorded SHA256. An archive digest never substitutes for a file digest.

Thirty-five Conda file digests match. Five remain unknown because there is
no directly usable installed-file digest: `libbz2.so.1.0.8`, `libgcc_s.so.1`,
`libgomp.so.1.0.0`, `libuuid.so.1.3.0` and `libz.so.1.2.13`. These are not
reported as mismatches or as package-byte verification successes.

For system paths, `dpkg-query -S PATH` identified `libc6:amd64`; `dpkg-query
-W` supplied version, architecture and installation status. All seven exact
file paths match their entries in the local package's `.md5sums` file. The
receipt preserves that manifest's identity and expected/observed MD5 values.
MD5 is only local consistency evidence, not cryptographic authenticity.
All file and referenced metadata identities were rechecked after collection.

## Remaining Work

Local metadata can be stale or modified and is not trusted upstream package
authentication. These observations do not establish a complete dependency
closure, verified acquisition recipe, compatible alternate runtime, source/
license attribution for statically linked components, security clearance or
redistribution rights. In particular, the OpenMP library observed in the
worker is supplied by the base Conda installation, outside the wheel lock.

The next packaging step must preserve this provenance while validating
separately acquired base-runtime artifacts and the resulting fresh runtime;
do not upgrade the running environment or silently replace historical bytes.
No installation, native process intervention, scientific result or benchmark
configuration changed during attribution. Job 22337 still requires terminal
accounting and independent full-workflow admission.
