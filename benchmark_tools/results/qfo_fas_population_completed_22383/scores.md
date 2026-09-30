# Retained FAS Eligible-Population Recount

Six unchanged recount rows are reused from timed-out attempt 22382;
FastOMA and OrthoMCL are freshly recounted. The timeout remains retained.
The continuation checks current input stability, not the original parser's
historical file identity or a missing original final stability pass.

Generated from a complete eight-method recount. Completion bounds assume
uncomputed eligible scores have hypothetical values in [0,1]. They are not
confidence intervals or replacement native FAS scores. Database hashes are
new retained identities, not independent historical checksum bindings.

| Method | Eligible pairs | Precomputed | Uncomputed | Precomputed mean | Hypothetical full-mean bound |
| --- | ---: | ---: | ---: | ---: | ---: |
| OrthoHMM high sensitivity | 9,032,719 | 6,902,660 | 2,130,059 | 0.790598 | 0.604163 - 0.839979 |
| OrthoHMM phylogeny satellite_v2 | 5,959,560 | 5,653,739 | 305,821 | 0.773667 | 0.733966 - 0.785282 |
| OrthoFinder 3.1.5 full | 14,215,382 | 11,819,171 | 2,396,211 | 0.728795 | 0.605946 - 0.774511 |
| OrthoFinder 3.1.5 sequence-only checkpoint | 163,277,439 | 29,303,161 | 133,974,278 | 0.691379 | 0.124081 - 0.944612 |
| SonicParanoid 2.0.9 | 15,248,739 | 12,675,597 | 2,573,142 | 0.739257 | 0.614511 - 0.783256 |
| ProteinOrtho 6.3.6 | 4,695,385 | 4,679,990 | 15,395 | 0.813883 | 0.811215 - 0.814493 |
| FastOMA 0.3.5 final orthologous groups | 15,008,180 | 7,616,298 | 7,391,882 | 0.775904 | 0.393753 - 0.886276 |
| OrthoMCL 1.4 | 14,378,089 | 13,127,208 | 1,250,881 | 0.742282 | 0.677704 - 0.764703 |
