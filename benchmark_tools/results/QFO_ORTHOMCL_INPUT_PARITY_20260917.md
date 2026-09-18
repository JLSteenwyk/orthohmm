# Retained OrthoMCL Input Parity

The retained QfO OrthoMCL 1.4 combined FASTA and genome map match the frozen
original OrthoHMM input set: 78 proteomes, 976,504 sequence IDs, and 976,504
identical case-sensitive sequence length/SHA-256 records. There are zero
sequence differences, missing or duplicate IDs, or incorrect species
assignments in `all.gg`. FASTA wrapping and header descriptions are not
sequence content; original first-token identifiers are preserved.

Evidence: `qfo_orthomcl_input_parity_20260917.json`, SHA-256
`d686d0ae252480e385ef2da5d274e8e9f84db60f776a5fe5721065dc1b757a0c`.
The report records hashes of all original FASTAs, combined FASTA, genome map,
frozen manifest and audit source. Inputs are checked before and after audit.

`run_orthomcl_qfo.slurm` names `work/all.fa` and `work/all.gg` as the combined
inputs. This audit checks those retained files, not the binary BLAST database
or an independently authenticated historical execution trace. It does not
erase the documented 53 sequence-specific BLAST failures. It demonstrates
parity with the original release, including its Xenopus incompatibility with
the scorer reference, not compatibility with the corrected release. No
sequences, predictions, scores or existing checkpoints were modified.

## Reproduction

Use a fresh output path; the command refuses to overwrite an existing report.

```sh
python benchmark_tools/audit_orthomcl_input_parity.py \
  --prepared benchmark_tools/results/qfo_factorial_prepared_20260917.json \
  --combined qfo_benchmark/results/orthomcl_1_4/work/all.fa \
  --genome-map qfo_benchmark/results/orthomcl_1_4/work/all.gg \
  --output /tmp/qfo_orthomcl_input_parity.json
```

Twelve focused tests cover exact and changed sequences (including case and
ambiguous-residue changes), reordered and wrapped records, incomplete,
unknown and duplicate combined IDs, malformed and incomplete genome maps,
wrong species assignments, duplicate original IDs and empty proteomes.
Combined with the existing OrthoFinder parity tests: 21 passed.

Remaining comparator input audits: SonicParanoid, Proteinortho and FastOMA.
Their converted pair files alone are not proof of native sequence parity.
