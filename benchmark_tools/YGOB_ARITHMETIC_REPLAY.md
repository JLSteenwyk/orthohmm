# Frozen YGOB Arithmetic Replay

This standalone workflow checks the frozen novel-taxon transfer results from
identifier-free sufficient statistics. It does not run OrthoHMM/OrthoFinder,
reconstruct the curated reference, re-admit native outputs, test family
independence or measure controlled resources. Scientific results and settings
are unchanged. It is not the complete publication release.

## Inputs and Environment

Use Python 3.12.3 and NumPy 2.2.6 for the recorded replay. The standalone script
imports no project modules and requires no Git, native orthology tools, raw
proteomes, YGOB pillar memberships or historical workstation files. NumPy is
the sole nonstandard Python import; compiled/OS dependencies are not bundled
or asserted hermetic. Use the accompanying requirements in a private environment,
not an upgrade of a shared scientific interpreter.

| Input | Bytes | SHA256 |
| --- | ---: | --- |
| `ygob_frozen_results_20260916.json` | 7,978 | `3927d81b1fd851eef3c8ce1343679558ef4f93bcfe2e643435e369587a44a4ef` |
| `ygob_sufficient_counts_20261002.json.gz` | 77,438 | `62885bb85eb357d97e55d053134d262e1ebb1e2b088a37df96986f98f8620d66` |

The compressed observation decodes to 590,099 bytes, SHA256
`2b931b0d11dba223f5f4ee3617eb88ae9a019f1fecb71809472d2ea812d6b89a`.
It preserves all 10,250 pillar sizes and four methods' 205,000 integer cells:
TP, twice the allocated FP, FN, covered genes and exact-recovery indicator.
The half-FP representation is exact, not rounding. One hash identifies the
ordered original pillar/size signature; the labels themselves are omitted.
No gene IDs, pillar names, sequences, reference memberships or predictions
are included. The summary retains historical file paths only as provenance.

These are project-computed numerical observations from the curated YGOB
homolog groups, not a replacement distribution of its raw data. Preserve YGOB
scientific attribution and the original study's acquisition/rights limitations;
the project license does not grant new upstream permissions. Identifier removal
does not itself establish redistribution clearance for an arbitrary future bundle.

## Execute

Keep the script, summary and compressed counts together after transfer. From
any directory, select a fresh output file whose parent exists:

```bash
python -I -B /copied/reproduce_ygob_validation.py reproduce \
  --snapshot /copied/ygob_sufficient_counts_20261002.json.gz \
  --snapshot-sha256 62885bb85eb357d97e55d053134d262e1ebb1e2b088a37df96986f98f8620d66 \
  --summary /copied/ygob_frozen_results_20260916.json \
  --output /fresh/ygob-arithmetic.json
```

Pin the transferred script/archive identity independently before executing it;
an internal index alone is not authenticity. The NumPy version is checked.
The script bounds decoded input and bootstrap batch memory, refuses overwrite,
preserves a failed reproduction report and checks input bytes before and after.
Its batch size defaults to 128 and may only decrease to a positive integer.
Outputs belong outside a byte-verified bundle so its inventory stays unchanged.

Point verification uses integer half-counts and checks all four methods' 12
ratios against the admitted summary and independent enumerated global counts.
Coverage/exact-recovery totals come from the retained pillar observations.
Prediction-group projection metadata is checked internally, not reconstructed
from omitted memberships. The sequence-only checkpoint remains diagnostic.

Interval verification regenerates the original 20,000 paired multinomial pillar
replicates, PCG64 seed 20260917, recomputes micro ratios in each replicate, and
checks both OrthoHMM-versus-full-OrthoFinder contrasts across F1, precision and
recall. All six effects and 24 nominal/Bonferroni bound values must agree within
1e-12, using linear quantiles and the original six-endpoint adjustment. It is
independent count arithmetic but the same NumPy RNG/quantiles, not a new
statistical engine or proof of pillar exchangeability.

The primary F1 difference remains -0.084870 percentage points, adjusted interval
[-0.622528, 0.445222]. This does not establish equivalence or superiority.
Novel-taxon transfer retains homolog-family overlap and cannot be promoted to
family-disjoint or unrestricted generalization by reproducing its arithmetic.

## Projection Provenance

The `export` command reads only the separately pinned original full result
(7,945,362 bytes, SHA256
`5104e000d0c3bb8687103646be8011fbc97f8b05f7eb0a925bbfdd1e76200e83`)
and the existing summary. It verifies their metadata, each method's original
ordered pillar/size signature and all global point/coverage identities before
writing the compact observation. It does not access source predictions or
reference memberships. That full result is not needed for standalone replay.

The first projection attempt failed before writing because the summary adds
the `full_results` provenance field absent from the full report. The exporter
now checks that exact one-field schema relationship and all remaining common
metadata, rather than silently ignoring arbitrary differences. No native
experiment, original result or endpoint was changed. Historical failures and
the first arithmetic replay remain separate evidence.
