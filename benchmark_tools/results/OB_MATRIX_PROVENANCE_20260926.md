# SonicParanoid and Proteinortho OrthoBench Provenance

The [native-table audit](ob_matrix_provenance_20260926.json) reproduces both
retained benchmark partitions, including evaluator-only singleton padding.
It does **not** pass all input-identity checks: the Proteinortho staged inputs
differ from the frozen benchmark in three files. The failed gate and original
receipt are preserved, not replaced by a passing normalization assumption.

| Method | Native groups | Native assigned genes | Added singleton groups | Native-to-benchmark partition |
|---|---:|---:|---:|---|
| SonicParanoid 2.0.9 | 14,849 | 206,466 | 44,912 | Exact |
| Proteinortho 6.3.6 | 23,588 | 184,897 | 66,481 | Exact |

Both normalized partitions cover 251,378 input gene IDs exactly once. Native
species columns, declared row counts, membership uniqueness and gene identity
were checked. Padding preserves pair co-membership; it does not increase native
assigned-gene coverage. The prediction hashes match previously admitted scores.

## Input Difference

All 12 SonicParanoid staged FASTAs match the frozen input hashes. Nine of twelve
Proteinortho staged FASTAs match; the remaining three retain identical gene
IDs/counts but have these sequence differences:

| Species | Changed proteins | Removed `*` symbols | Changed genes in full RefOG memberships |
|---|---:|---:|---:|
| Danio rerio | 6 | 21 | 0 |
| Homo sapiens | 90 | 322 | 2 |
| Mus musculus | 81 | 526 | 0 |
| Total | 177 | 869 | 2 |

Every changed sequence is exactly explained by deleting all `*` characters,
including internal occurrences, from its frozen counterpart. This is an
observed transformation, not attribution of who changed the files or when.
The changed reference genes are `ENSP00000486295` and `ENSP00000487059`, both
in `RefOG023.txt`. Exposure uses full reference membership before low-certainty
exclusions. [Per-protein differences and sequence hashes](ob_proteinortho_input_differences_20260926.json).

This does not quantify the accuracy effect: changed proteins outside RefOGs
can affect false positives and clustering, and unchanged IDs do not establish
unchanged inference. Do not describe Proteinortho's historical OrthoBench row
as a verified byte-identical-input comparison. Its score remains the retained
development-exposed result, with this additional comparability limitation.

## Execution Evidence and Gaps

SonicParanoid's native run-info records version 2.0.9, 32 threads, DIAMOND
very-sensitive mode, MCL inflation 1.5, and the expected input/output paths.
Its log reports 1219.933 seconds total elapsed time. This is tool-reported,
not an independent end-to-end timing record; no CPU time or peak RSS is claimed.
Exact invocation argv was not recovered from these records.

Proteinortho's info file retains two identical command invocations on April 20,
with `-cpus=32 -project=orthobench`, version 6.3.6 and parameter vectors.
The main log reports completion and DIAMOND 2.1.12. The multiple invocations
are retained rather than presented as one uninterrupted run. Per-invocation
outcomes, elapsed/CPU time and peak RSS remain unverified. Staged-file evidence
does not independently attest historical consumption or binary identity.

## Reproduction

```bash
python -m benchmark_tools.audit_ob_matrix_provenance \
  --root . --output /tmp/ob-matrix-provenance-new.json
# Expected exit status 1: Proteinortho byte-identity gate fails; JSON is retained.
python -m benchmark_tools.explain_ob_input_differences \
  --report /tmp/ob-matrix-provenance-new.json \
  --reference benchmark_tools/results/orthobench_paired_uncertainty_20260916.json \
  --output /tmp/ob-input-differences-new.json
```

Twenty-six focused input, matrix, metadata and difference tests pass. No native
inference was rerun, scores changed, input files repaired, or controlled timing
claim made. Broader provenance and publication requirements remain open.
