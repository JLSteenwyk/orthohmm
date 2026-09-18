# Original QfO Comparator Input Parity

This audit checks the original release only. It does not resolve the Xenopus
release mismatch with the scorer reference, establish corrected-release
accuracy, or justify replacing the inputs of existing runs.

## Retained Evidence

| Method | Retained representation audited | Result | Boundary |
| --- | --- | --- | --- |
| OrthoFinder full | Internal Species FASTAs plus species/sequence maps | All 976,504 sequences, 78 species, exact sequence parity | Not an authenticated historical execution trace |
| OrthoFinder sequence-only | Same historical run's MCL checkpoint | No separate input audit here | Checkpoint derivation/conversion needs its own evidence |
| OrthoMCL 1.4 | Combined `all.fa` and `all.gg` | All 976,504 sequences identical; complete correct species map | Not binary BLAST database validation; 53 failed queries remain |
| SonicParanoid 2.0.9 | `results/sonicparanoid/run/input/*.fasta` | All 78 complete file hashes identical | Staged copies, not downstream native representations |
| Proteinortho 6.3.6 | `results/proteinortho/run/input/*.fasta` | All 78 complete file hashes identical | Staged copies, not binary search database validation |
| FastOMA 0.3.5 | `results/fastoma/run/input/proteome/*.fa` | All 78 complete file hashes identical | Staged copies; UniProt ID transformation not audited here |

Paths in the table are relative to `qfo_benchmark`. The reference for all
comparisons is the original frozen input manifest
`qfo_factorial_prepared_20260917.json`, SHA-256
`706b07c91e9a130dae229837641a7ad62d7d36a09679e5e0daa0959e182b7d64`.
This original sequence universe is also the frozen OrthoHMM ablation input.
The full-file matches include headers and wrapping, not just gene totals.
All 234 staged files are regular copies with distinct resolved paths from
the original inputs; none is reported as an original-input symlink.

## New Three-Tool Audit

Evidence `qfo_comparator_staged_input_parity_20260917.json`, SHA-256
`e5f96229368e651967dd2f0e4b0000e009299e33fc2ee5092d0d5cafea86e56f`.
The audit records each original/staged file identity, checks expected names,
rejects extra common-suffix FASTAs, checks files before/after comparison and
records symlink/shared-path status. No data or predictions were modified.

SonicParanoid's retained `run.log` explicitly names the audited staging
directory and 78 proteomes. FastOMA's retained `run.log` names its audited
input folder, version 0.3.5 and `--fasta_header_id_transformer UniProt`.
The FastOMA wrapper copies `.fasta` to `.fa` without rewriting content.
Proteinortho's log lists checks of input FASTA basenames. These logs support
the selected paths but do not cryptographically bind historical execution
to current file contents. FastOMA's supplied OrthoFinder species tree is a
separate methodological limitation, not removed by sequence parity.

Earlier detailed reports:
`QFO_ORTHOFINDER_INPUT_PARITY_20260917.md` and
`QFO_ORTHOMCL_INPUT_PARITY_20260917.md`.

## Reproduction

```sh
python benchmark_tools/audit_qfo_staged_inputs.py \
  --root . \
  --prepared benchmark_tools/results/qfo_factorial_prepared_20260917.json \
  --output /tmp/qfo_comparator_staged_inputs.json
```

Use a fresh output filename. Eight focused tests cover identical and changed
bytes, explicit symlink reporting, extra/missing FASTAs, duplicate expected
names and original-hash mutation. The combined staged-input, OrthoMCL and
OrthoFinder parity suites pass 29 tests.

Next: complete corrected-archive mapping and native sequence compatibility
checks before freezing any separately named corrected-input experiment.
Keep current original-release scores labeled as release-limited. Agreement
of retained input contents does not prove equal effects of the release
mismatch on different methods or invariant rankings after correction.
