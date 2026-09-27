# Matched-Graph Stage Trace

Post-hoc descriptive diagnostic; no causal or confirmatory inference.

All values below are equal-dataset mean overlap F1 percentages.
Direct hits, graph edges and component closure are not native ortholog predictions.

| Condition | Stage | HMM | DIAMOND |
|---|---|---:|---:|
| baseline | direct_hits | 94.6336 | 95.0810 |
| baseline | rbnh_edges | 94.6287 | 94.8262 |
| baseline | rbnh_components | 95.5616 | 95.3749 |
| baseline | initial | 95.7240 | 95.3749 |
| baseline | multipass_edges | 94.6287 | 94.8262 |
| baseline | multipass_components | 95.5616 | 95.3749 |
| baseline | multipass | 95.7240 | 95.3749 |
| baseline | final | 95.7240 | 95.3749 |
| divergent | direct_hits | 46.0226 | 45.4181 |
| divergent | rbnh_edges | 45.8713 | 43.9013 |
| divergent | rbnh_components | 64.2299 | 49.5592 |
| divergent | initial | 59.0573 | 49.3443 |
| divergent | multipass_edges | 45.8746 | 43.9285 |
| divergent | multipass_components | 64.2865 | 49.6710 |
| divergent | multipass | 59.0678 | 49.4565 |
| divergent | final | 59.0789 | 49.4565 |
| divergent_turnover | direct_hits | 45.9566 | 44.6396 |
| divergent_turnover | rbnh_edges | 45.4281 | 42.4456 |
| divergent_turnover | rbnh_components | 65.6061 | 49.1586 |
| divergent_turnover | initial | 59.9822 | 48.8569 |
| divergent_turnover | multipass_edges | 45.4744 | 42.5059 |
| divergent_turnover | multipass_components | 65.6586 | 49.3164 |
| divergent_turnover | multipass | 60.1303 | 49.0173 |
| divergent_turnover | final | 59.5801 | 48.8969 |
| missing20 | direct_hits | 94.7533 | 95.0905 |
| missing20 | rbnh_edges | 94.6564 | 94.8466 |
| missing20 | rbnh_components | 95.4935 | 95.3321 |
| missing20 | initial | 95.6338 | 95.3321 |
| missing20 | multipass_edges | 94.6564 | 94.8466 |
| missing20 | multipass_components | 95.4935 | 95.3321 |
| missing20 | multipass | 95.6338 | 95.3321 |
| missing20 | final | 95.6338 | 95.3321 |
| taxon_count_control | direct_hits | 94.7516 | 95.0499 |
| taxon_count_control | rbnh_edges | 94.4312 | 94.4818 |
| taxon_count_control | rbnh_components | 95.3561 | 95.1231 |
| taxon_count_control | initial | 95.7460 | 95.1231 |
| taxon_count_control | multipass_edges | 94.4312 | 94.4923 |
| taxon_count_control | multipass_components | 95.3561 | 95.1637 |
| taxon_count_control | multipass | 95.7460 | 95.1637 |
| taxon_count_control | final | 95.5191 | 95.1637 |
| turnover | direct_hits | 84.9328 | 85.4528 |
| turnover | rbnh_edges | 87.2062 | 87.2994 |
| turnover | rbnh_components | 85.8125 | 85.7799 |
| turnover | initial | 86.6471 | 85.9741 |
| turnover | multipass_edges | 87.2062 | 87.2994 |
| turnover | multipass_components | 85.8125 | 85.7799 |
| turnover | multipass | 86.6471 | 85.9741 |
| turnover | final | 84.9816 | 84.5499 |
| uneven_taxa | direct_hits | 93.7331 | 94.5599 |
| uneven_taxa | rbnh_edges | 93.4165 | 93.9382 |
| uneven_taxa | rbnh_components | 94.7504 | 94.8518 |
| uneven_taxa | initial | 94.7504 | 94.8518 |
| uneven_taxa | multipass_edges | 93.4165 | 93.9382 |
| uneven_taxa | multipass_components | 94.7504 | 94.8518 |
| uneven_taxa | multipass | 94.7504 | 94.8518 |
| uneven_taxa | final | 94.7504 | 94.8518 |
| overall | direct_hits | 79.2548 | 79.3274 |
| overall | rbnh_edges | 79.3769 | 78.8199 |
| overall | rbnh_components | 85.2586 | 80.7399 |
| overall | initial | 83.9344 | 80.6939 |
| overall | multipass_edges | 79.3840 | 78.8339 |
| overall | multipass_components | 85.2742 | 80.7843 |
| overall | multipass | 83.9571 | 80.7386 |
| overall | final | 83.6097 | 80.5180 |

## Limits

- Post-hoc, development-exposed diagnostic, not a new confirmatory comparison.
- Overlap scores for search hits, edges and component closure are not ortholog prediction accuracy.
- A pair without direct search evidence can be recovered through indirect graph paths.
- Stage localization does not isolate hit identity, score, ranking or clustering effects causally.
- Five seeds per condition share histories across conditions; no new intervals or tests.
