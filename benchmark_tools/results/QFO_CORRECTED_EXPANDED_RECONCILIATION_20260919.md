# Corrected Expanded-Candidate Reconciliation

Array task21760_1 (raw job21862) completed0:0 in1:54:53 on the shared
host; independent native admission21762 completed0:0 in2:22. The cell is
p0_c1_r1: initial HMM search retained, multi-sequence profile expansion off,
candidate expansion on, inferred phylogenetic reconciliation on. It is not
the full publication p1_c1_r1 configuration.

The [unchanged native admission receipt](qfo_corrected_factorial_native_admission_21762.json)
has SHA-256 `27d359b059d3c7589e8f85adf03a00fc5e53d57e6b80f4e8b75727f0749841a8`.
The admission checks146013recorded artifacts and partition/native-pair
integrity. This follow-up additionally rehashed23distinct referenced files,
including execution records, native pairs, metrics, tree and helper sources.
It does not claim a new independent reconstruction of every tree/event.

All984137genes are preserved across368209root HOGs from353638candidate
families.7358candidate families split and no cross-source merges are
reported. The membership policy supports29089of40690constraints and detaches
11601, involving26477genes; it reports185610removed ortholog pairs and
11443added root HOGs. These are algorithmic decisions, not truth labels.

Conversion21768 completed0:0 in3:29. Its fresh independent native recheck
matches the admission receipt exactly. All5977100native phylogenetic pairs
remain after reference mapping, with zero mapping loss and no root-HOG
clique expansion. [Conversion receipt](qfo_corrected_factorial_pairs_21768.json)
SHA-256: `10bac6240a02d188faaf3f2a6a7a89e9b3080274d6ab110afcb1043d9904907a`.
Both91,974,760-byte pair files rehash to
`98726096cba488eb24d64d03fde77a73fcd51653a216e8ae347fb39277c22e13`.
Scoring21779 is running; independent assessment21780 remains dependent.
No new accuracy value is admitted here.

## Exploratory Mechanism Check

[Generated native-output table](qfo_corrected_p0_reconciliation_comparison_20260919.tsv)
and [machine-readable comparison](qfo_corrected_p0_reconciliation_comparison_20260919.json)
compare p0_c0_r1 with p0_c1_r1. The JSON SHA-256 is
`e6dd4f273ad89012891aad83dd3fb1a206b6d625634b0f49c84a91b7b8bd7ee2`.
Candidate counts decrease394328to353638, root HOGs397041to368209, and
native pair counts increase5113820to5977100. No accuracy conclusion follows
from these volumes.

The inferred species trees also differ. The runs use18versus25families
for species-tree estimation. Both trees contain the same78species and
76nontrivial rooted clades. They share38clades, with38unique to each, giving
rooted symmetric difference76. DendroPy5.0.8's comparison agrees exactly
with an independent set comparison of descendant-leaf clades; branch lengths
are ignored. Neither tree is ground truth.

This is a post-outcome exploratory diagnostic, not an added primary endpoint
or basis for tuning. The candidate-expansion contrast includes downstream
species-tree estimation as part of its end-to-end effect. It cannot be
interpreted as candidate expansion holding the realized tree constant, or
as evidence that tree changes caused a particular accuracy change. Such a
mechanistic claim would need separately specified fixed-tree controls and
appropriate truth; no such extra run is claimed here.

## Reproduction

```bash
/home/bizon/anaconda3/bin/python benchmark_tools/compare_corrected_reconciliation.py \
  --left benchmark_tools/results/qfo_corrected_factorial_native_admission_21761.json 97341e1b9ef6ac36e6b8c329aa3b5161690ccf00e617895f1d4ccb6ea127b16a \
  --right benchmark_tools/results/qfo_corrected_factorial_native_admission_21762.json 27d359b059d3c7589e8f85adf03a00fc5e53d57e6b80f4e8b75727f0749841a8 \
  --output FRESH_COMPARISON.json
```

Eleven focused tests cover topology/branch-length invariance, changed and
unresolved topologies, unequal taxa, checksum/provenance/count mismatches,
changed tree bytes and non-admitted inputs. The final evaluator exactly
reproduces the retained JSON. Shared-host elapsed times are incremental
reconciliation/conversion stages, not dedicated end-to-end comparisons.
