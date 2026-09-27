# Full OrthoBench Scoring Outside The Checkout

The [portable scoring receipt](portable_ob_score_20260927.json) connects the
verified upstream acquisition, exported frozen scientific readers and patched
reader environment to the full retained OrthoBench predictions. It recomputes
scores from gene memberships and raw reference files, not retained counts.

All 251,378 genes are covered exactly once in 59,770 root groups. The complete
score object, including all 70 reference-family records, matches the admitted
full recovery run exactly:

| Metric | Reproduced value (%) |
| --- | ---: |
| Weighted RefOG F1 | 74.10607351873405 |
| Precision | 81.77045380181866 |
| Recall | 67.75533630827974 |

There are 15 exact RefOGs. Historical predictions and results are unchanged.
This is post-inference scoring, not another native inference run or new
evidence of biological accuracy.

## Execution And Checks

The controller revalidates the separately acquired upstream checkout at
`872d6f30592ab5ff837224db16a514b3f2bb916a` against Git blobs and the frozen
input manifests. It copies 12 FASTA files, 70 reference files, 11 low-certainty
files and the admitted root-group TSV to a fresh private directory. The
31-module reader export is also copied without changes. The worker checks
the complete 94-file data inventory, hashes, all FASTA IDs, reference IDs and
strict root-group membership before invoking the unchanged frozen scorer.

Successful root: `/tmp/orthohmm-portable-ob-score-20260927-v2`. The worker used
`/tmp/orthohmm-reader-patched-20260927/venv/bin/python -I -B` with a clean
environment. Its file trace contains neither the original repository prefix
nor the main installation's site-packages prefix. The worker does not read
the expected score; the controller compares results after worker completion.
Input and exported-source identities are checked again afterward.

The initial attempt stopped before scoring because a new validation check
treated an empty line in a low-certainty reference as an unknown gene ID.
Diagnostic inspection confirmed 251,378 FASTA IDs, no missing nonempty
reference IDs, and only the empty-string discrepancy. The corrected check
ignores the empty string only for ID-universe validation; original reference
contents and frozen scorer inputs are unchanged. The failed directory,
`/tmp/orthohmm-portable-ob-score-20260927`, retains its log and failure receipt.
The second attempt used a new directory; no automatic retry or overwrite.

Forty-seven focused tests pass across the portable runner, acquisition
verifier, frozen scorer and strict root-group parser. Nine receipt identities
and all 94 staged data files were independently rechecked after execution.

## Reproduction

After the separately documented acquisition, reader export and reader-only
installation steps, run from the project checkout:

```bash
python -m benchmark_tools.reproduce_portable_ob_score \
  --repo . \
  --acquisition /acquired/Open_Orthobench \
  --readers /exported/independent-readers \
  --python /patched/reader/venv/bin/python \
  --output /fresh/portable-ob-score \
  --report /fresh/portable-ob-score.json
```

The controller still requires the admitted historical run at its documented
repository-relative location. The copied worker uses only its relocated
inputs and exported modules; it can be invoked separately with `--worker`
and a fresh `--output`. This distinction is intentional: full release
restoration and a combined acquisition-to-inference workflow remain open.
Raw benchmark data are private acquisition copies, not redistribution-cleared
release assets. No new official upstream scorer execution, timing comparison,
cross-host validation, QfO scoring or publication-readiness claim is made.
