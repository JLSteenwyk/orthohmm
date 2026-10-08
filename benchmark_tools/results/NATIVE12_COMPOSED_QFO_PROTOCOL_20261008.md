# Final Native QfO Conversion

This prospective workflow consumes the explicit `native12_composed_terminal_review_v1`
for the frozen QfO identity 12, `p1_c1_r1`, inference job 24036. It does not
translate that result to an original terminal-review schema. The inference,
controller, terminal reviewer and their bound sources remain unchanged.

`native12_composed_review_binding.py` requires the actual inference-bound
request digest, new review identity, successful native and reviewer accounting,
reviewer source/batch already bound by that request, and one retained owned held
submission plus one release. The future reviewer receipts must use
`native12_composed_review_submission_v1` and `native12_composed_review_release_v1`,
with the original request/worker/batch references and actual raw held envelope.
The reviewer stdout must contain the exact returned review file record.
The full resource/environment/runtime and independently validated scientific
output components and their retained evidence are rehashed before conversion.
Historical OS-inventory inequality and retained earlier review failures are
explicit; no continuous integrity, isolated timing or original-schema claim is made.

`prepare_native12_composed_qfo_pairs.py` requires a scheduled two-CPU, 32 GiB
conversion and a fresh separate destination. It uses the unchanged native
conversion, mapping filter and coverage kernels. Its input is specifically
`native/orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv`, already present in
the scientific review's checked files. The row count must match that review's
`phylogeny.native_pair_rows`. Group cliques and pre-clustering graph edges are
not substitutes. Accession normalization must be injective across the complete
input universe; zero mapping loss is required. Coverage retains every input
accession in its denominator, including singletons.

The batch accepts request path/digest, review path/digest, held path/digest,
release path/digest, actual review job ID, and converter source digest, in that
order. It has a 12-hour scheduler limit and no requeue. Submit only after actual
successful terminal review; retain/check the owned held envelope, commit it,
then capacity-check and release once. No production conversion is currently
submitted. A failed conversion retains preflight, partial predictions and its
failure report; no silent retry or imputation is authorized.

Conversion produces `native12_composed_qfo_conversion_v1`, status
`native12_composed_qfo_pairs_prepared_unscored`, participant
`ohmm_qfo_full_native_p1_c1_r1`. It is not accuracy admission. Prospective
six-endpoint assessment and independent admission still need their compatible
consumers; use the unchanged frozen QfO 2020 and FAS kernels, no endpoint changes
or resume. Failed inference/review means no converted prediction or score.

Focused tests use real small native-pair/filter/coverage kernels, explicit
synthetic scheduler/review contexts, schema and envelope mutation tests,
producer-stdout binding, empty predictions, full denominator, unmapped
singletons, native-count mismatch, clique rejection, and retained postflight
failure. These tests are not a production review or scientific score.
