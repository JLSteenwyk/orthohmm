# Selected Wheel Native Dependency Inventory

The [machine-readable inventory](integrated_wheel_elf_20260927.json) inspects
the same 12 exact wheel artifacts as the integrated workflow notice export.
Both installation-bound notice inventories were recomputed before inspection;
wheel identities were checked before and after each scan. Every regular
member was checked for ELF magic, rather than relying on filename suffixes.
Sixty ELF objects were found, matching all 60 earlier filename candidates.
No extra ELF object or non-ELF filename candidate was found in this set.

GNU [readelf](https://sourceware.org/binutils/docs/binutils/readelf.html)
`--wide --dynamic` reads each object from an isolated temporary file without
loading or executing it. The receipt pins the inspector binary/version,
source, wheel/member digests and raw dynamic-section text. Diagnostics stop
the scan. Duplicate, unsafe and nonregular archive members are rejected;
archive paths never determine temporary extraction paths.

## Observations

- 182 NEEDED edges were recorded. Eleven have basename/SONAME candidates
  among the selected wheels, covering eight distinct dependency names.
- Bundled candidate names include igraph, Leiden, OpenBLAS, libgfortran,
  libquadmath, libgomp, libxml2 and liblzma artifacts. These are filename
  observations, not validated component versions, license or security findings.
- Twelve distinct names lack selected-wheel candidates: `ld-linux-x86-64.so.2`,
  `libc.so.6`, `libdl.so.2`, `libgcc_s.so.1`, `libgomp.so.1`,
  `libgomp.so.1.0.0`, `libm.so.6`, `libpthread.so.0`, `librt.so.1`,
  `libstdc++.so.6`, `libtbb.so.12` and `libz.so.1`.
- Seven objects declare RPATH entries relative to `$ORIGIN`; no RUNPATH
  entry was observed. The report preserves the exact strings.
- The TBB declaration belongs to `numba/np/ufunc/tbbpool`'s extension.
  This does not establish that the workflow loads that backend or requires
  a new TBB installation. No package or runtime was modified.

## Scope And Next Steps

This is not dynamic-loader resolution: candidate matching ignores environment
separation, load scope and search order. Missing selected-wheel candidates
may be supplied by the operating system or belong to unused optional
extensions. Static linking and `dlopen` dependencies remain outside this
scan, as do external MAFFT/FastTree executables and Python/OS libraries.
The notice lists are retained alongside objects, but are not a completed
component-to-license mapping. Inspect loaded components and their source/
notice obligations before making runtime portability or redistribution claims.

The subsequent [graph source acquisition](GRAPH_SOURCE_CANDIDATES_20260928.md)
retrieves hash-verified igraph/leidenalg source candidates and identifies
external library provenance still missing; it does not close this clearance gap.

Forty-one focused ELF/notice tests pass, including content-based detection,
temporary-file cleanup, changed-wheel rejection, unsafe/duplicate/nonregular
members, malformed tags, inspector diagnostics and unresolved candidates.
The earlier full-suite result predates this new scanner and remains separate.

```bash
python -m benchmark_tools.inventory_wheel_elf \
  --inventory benchmark_tools/results/integrated_inference_notices_20260927.json \
  --inventory benchmark_tools/results/integrated_reader_notices_20260927.json \
  --output /fresh/integrated-wheel-elf.json
```

No scientific score, running job, historical artifact, legal clearance or
publication-readiness status was changed.
