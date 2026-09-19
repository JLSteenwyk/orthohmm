# Corrected Profile-Refined Reconciliation

Array task21760_2 (raw job21888) completed0:0 in1:14:07 on the
shared host. Independent native admission21763 completed0:0 in2:20.
The cell is p1_c0_r1: HMM profile refinement on, candidate expansion
off, inferred phylogenetic reconciliation on. This is not the full
publication p1_c1_r1 configuration.

The [unchanged admission receipt](qfo_corrected_factorial_native_admission_21763.json)
has SHA-256 `36f3c8a3b9cebf8f474b3cb956d3e3e79cebb5243a4baa52e89a1d0a5c23d317`.
It records157209artifacts and validates native partition/pair integrity.
This follow-up independently rehashed all23distinct path/size/SHA-256
records referenced in the receipt, including execution status, native
pairs, metrics, species tree and helper sources. This is not an independent
reconstruction of gene trees or evolutionary truth.

All984137genes are preserved:391908candidate families yield394559root
HOGs, with669source families split and zero cross-source merges. The native
phylogenetic output contains5113180pairs. These are structural counts, not
accuracy estimates; pairs must not be replaced by root-HOG cliques.

Conversion21770 completed0:0 in3:21. Its fresh native-admission recheck
matches the original receipt byte-for-byte. All5113180native pairs survive
reference mapping, with zero mapping losses. Both converted files contain
78,540,928bytes and have SHA-256
`45e3bd4ae7a35a517b6de4a599b2a336566862a5e915463a9c14892bdf907c61`.
The [unchanged conversion receipt](qfo_corrected_factorial_pairs_21770.json)
has SHA-256 `92ac5604612fbc52a9a5ffeea95f8ad6255d8652f480a0faebec098f33c0ac00`.
Independently rehashed all106distinct referenced file records, including
inputs, mapping, both converted pair files and the repeated admission.

Scoring21783 is running; independent score admission21784 remains
dependent. No accuracy result is admitted yet. The final reconciliation
task21760_3 is running.

Elapsed times above describe incremental stages on the shared host, not
dedicated end-to-end efficiency comparisons. The corrected eight-cell
uncertainty analysis21894 still awaits its complete input set. No default,
endpoint, bootstrap protocol or running executor was changed.
