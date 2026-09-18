# Retained OrthoFinder Input Sequence Parity

Audited the retained internal input FASTAs of the QfO full OrthoFinder3.1.5
run under `orthofinder_v3_diamond`, the inference directory recorded by
`qfo_orthofinder_v3_1_5_20260828.json` and used in the retained comparator
table. The selected working directory is
`run/input/OrthoFinder/Results_Jul15/WorkingDirectory` within that run.

## Result

All78proteomes and976,504sequences match the frozen original QfO inputs
exactly after restoring internal numeric IDs through `SequenceIDs.txt`.
No sequence-byte differences, missing input records, duplicated IDs or
species-assignment mismatches were found. `SpeciesIDs.txt` matches the
original78-file inventory. All expected internal `SpeciesN.fa` files exist
with no extras; every mapped record appears exactly once in its own species.

Machine evidence `qfo_orthofinder_internal_input_parity_20260917.json`
SHA-256:
`9ad7b0b32b240827d202338cfafdf56fc1086d15e161caff67b9b5f9095b3fea`.
It retains file hashes for the frozen input manifest, both ID maps, original
FASTAs and internal FASTAs, with per-species counts. These files were checked
before and after analysis. This is full sequence-content comparison, not
merely matching gene totals or similar downstream scores.

## Boundaries

The retained full-pipeline internal input set therefore contains the same
original Xenopus sequence and mapping limitations as the frozen original
OrthoHMM inputs. This is evidence against an input-release advantage in
this retained OrthoFinder input set. It does not prove the historical files
were never changed after execution, authenticate every inference/scoring
step, establish other tools' input parity, or quantify corrected-release
accuracy effects. Those claims require their own provenance and rerun checks.

The sequence-only comparator was derived from an MCL checkpoint; this audit
does not independently revalidate that conversion or transfer full-pipeline
semantics to the checkpoint. No original or corrected dataset is substituted.

## Reproduction

```sh
python benchmark_tools/audit_orthofinder_input_parity.py \
  --prepared benchmark_tools/results/qfo_factorial_prepared_20260917.json \
  --working qfo_benchmark/results/orthofinder_v3_diamond/run/input/OrthoFinder/Results_Jul15/WorkingDirectory \
  --output benchmark_tools/results/qfo_orthofinder_internal_input_parity_20260917.json
```

Nine tests cover mapping structure, descriptions containing colons, exact
and differing sequences, duplicate and unknown IDs, incomplete coverage,
incorrect species IDs and unsafe species filenames. No inference, scoring,
sequence normalization or input modification is performed.
