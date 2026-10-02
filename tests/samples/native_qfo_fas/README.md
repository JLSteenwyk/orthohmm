# Retained Native FAS Source

`fas_benchmark.py` and `LICENSE` are unchanged files from QfO's
`benchmark-webservice`, commit `c0854a96c1a0fd7f2a891d971af0863002fabc90`.
They remain subject to the upstream Mozilla Public License 2.0, not the
OrthoHMM project license. Full upstream notices, including its separate Darwin
system notice, are retained in `LICENSE`; no Darwin binary is included.

The [provenance inventory](provenance.json) pins exact bytes and original
location. Both files match public pinned Git blobs and the retained upstream
checkout; the scorer also matches the historical canonical source archive.
These are test references, not a replacement benchmark implementation.

Tests select these files explicitly and verify their fixed identities before
AST comparison. Only the exact native loader function is executed on small
synthetic lookup entries; module imports, FAS inference, annotation files,
databases, benchmark labels and sampled scores are not substituted or bundled.
The native query literal is compared with the auditor's query separately.

Production source-admission rules and historical source paths/receipts are
unchanged. This snapshot is not the full QfO runtime, raw scoring reproduction,
data-redistribution clearance or independent biological validation.
