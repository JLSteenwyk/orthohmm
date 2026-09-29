# Fixed-Release Proteome Reacquisition

All eleven fixed-release Ensembl proteomes used in the retained Three Kingdoms
input audit were freshly downloaded over HTTPS. Each matches the retained
compressed file exactly in byte count and SHA256. This adds independent
provider retrieval to the earlier compressed/raw/staged local consistency audit.

The [receipt](three_kingdoms_provider_reacquisition_20260928.json) records every
URL, response header, retrieval time, provider checksum text, local and fresh
file identity, auditor and checksum executable. All eleven files also match
their release directory's `CHECKSUMS` entry using GNU `sum -r`. That checksum
is a weak 16-bit BSD sum with rounded 1,024-byte block counts, not a substitute
for the separate full-download SHA256 comparison.

| Source | Codes | Exact fresh matches |
|---|---|---:|
| Ensembl release 100 | Anol, Dani | 2 |
| Ensembl Metazoa 47 | Amph | 1 |
| Ensembl Plants 47 | Arab, Sela, Zeam | 3 |
| Ensembl Plants 58 | Oryz | 1 |
| Ensembl Fungi 47 | Asfu, Scer, Scom, Spom | 4 |
| Moving UniProt URL | Xeno | Not attempted; historical release unresolved |

The initial checksum-only pass is retained locally. The second, expanded pass
adds exactly one fresh download per fixed source, with no automatic retry.
Fresh files are separate from historical inputs under
`benchmarks/work/three_kingdoms_provider_checksums_20260928_v2/`.
Thirty-one focused tests pass, including malformed checksum rejection, moving
source exclusion, exact/mismatched/short download handling, size limits,
retained partial bytes and refusal to overwrite existing outputs.

## Reproduce

With the retained source audit and compressed inputs available, this command
creates an unused output directory and verifies provider data against them:

```bash
python -m benchmark_tools.audit_three_kingdoms_provider_checksums \
  --source benchmark_tools/results/three_kingdoms_sources_20260918.json \
  --source-sha256 a79bf83e1ea28a9790e4597ea50a44d6409f4f36404a35e660d56a0b754fc1f3 \
  --output-directory /fresh/path/three-kingdoms-sources \
  --download-proteomes
```

Omit `--download-proteomes` for the weaker provider-checksum-only audit. The
source audit retains historical absolute paths; relocate its referenced data
explicitly rather than pretending this command is a standalone portable bundle.
Unresolved downloads and mismatches are distinct from matches in the report.

## Limits

Xenopus's moving `current_release` URL is deliberately not used to replace
historical bytes or infer their release date. Fresh provider equality does not
prove what every historical method consumed; method-specific input traces are
separate. No FASTA, staged sequence, BUSCO call, reference group or benchmark
score changed. In particular, the previously documented zebrafish stop-marker
removal is not reversed by this audit. Exact retrieval is not redistribution
permission: proteomes remain source-acquisition-only for new release packaging,
and BUSCO/derived-reference rights require separate review.
