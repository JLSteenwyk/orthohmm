# Relocatable Source Notice Supplement

The four source texts needed for the newly attributed bundled components now
have an executable export and offline verification workflow. This supplements,
without modifying, the existing 80 wheel-embedded notice candidates.

The supplement includes libxml2 `Copyright`, XZ `COPYING`, and GCC `COPYING3`
and `COPYING.RUNTIME`. It does not export the inspected libgomp header, any
binary, a source archive, distribution patch, or dataset. The four notice texts
total 42,535 bytes. Their 5,471-byte `SOURCE_NOTICE_INDEX.json` has SHA-256:

`055771521faae8a936fe0d96c18d0f46ecd4c11bb3186ccf729cd6ba05509432`

## Export And Relocation

`benchmark_tools.export_bundled_source_notices` reads each selected text directly
from the hash-bound nested source archive, checks it against the retained source
inventory, and checks the source RPM, archive-inventory and signature-receipt
bytes before finalizing the index. It does not rerun RPM signature verification
or infer that a retained receipt constitutes legal clearance. Missing or duplicate
texts, nonregular selected archive entries and changed evidence fail. A partial
failed directory is not overwritten or finalized as a valid export.

The [export receipt](igraph_source_notice_export_20260928.json) and
[relocation receipt](igraph_source_notice_relocation_20260928.json) bind the
actual output. All four payloads and the identical index verify at the relocated
root. Offline verification uses an externally supplied expected index digest;
it does not need the original source archive paths to exist. Additional files,
changed texts, symlinks and changed index bytes are rejected.

```sh
python -B -m benchmark_tools.export_bundled_source_notices \
  --inventory benchmark_tools/results/igraph_bundled_source_material_20260928.json \
  --output /fresh/source-notices \
  --receipt /fresh/source-notices-receipt.json
```

The Python `verify(directory, index_sha256)` function validates a relocated
copy against the separately retained export receipt. Production exports remain
local at `benchmarks/work/igraph_source_notice_supplement_20260928` and
`benchmarks/work/igraph_source_notice_relocated_20260928`.

43 focused source-export, wheel-export and notice-inventory tests pass. They
include verification after deleting synthetic original inputs, bad archive
entries, missing/extra export content, source corruption and changed manifest
bytes. This demonstrates packaging integrity, not license compatibility or
complete transitive coverage. Native build correspondence and final release
review remain open. No frozen runtime or scientific result changed.
