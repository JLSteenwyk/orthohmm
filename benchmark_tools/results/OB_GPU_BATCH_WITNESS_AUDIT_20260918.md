# OrthoBench GPU Batch Witness Audit

## Question

Could the [all-long-target routing defect](SEARCH_GPU_ROUTING_FIX_20260918.md)
explain the retained OrthoBench search results? An accepted hit whose target
has at most1998residues demonstrates a GPU-eligible candidate in its search
batch. Under one complete batch per directed species pair, one such hit
rules out the all-long condition for that direction. Lack of such evidence
would remain unresolved, not prove that the bug occurred.

## Result

The [machine-readable audit](ob_gpu_batch_witnesses_20260918.json) inspected
18,235,373retained hits across251,378proteins and12species:

| Property | Count |
| --- | ---: |
| Directed species pairs, including within-species | 144 |
| Directions with an eligible accepted target | 144 |
| Directions lacking an eligible accepted target | 0 |
| Retained hits to eligible targets | 18,026,537 |
| Fewest eligible-target hits in any direction | 18,635 |

For every direction, the report retains a deterministic example chosen by
target length, target ID and query ID. These example target lengths range
from15to51residues. No reference labels, accuracy scores or final groups
were used to select them.

The input cache SHA-256 was checked before trusted deserialization. All
12FASTAs were checked against the admitted family-trace inventory, and the
complete gene universe, species ownership and lengths were compared with
the cache. Source/input records were checked again after the scan. All7unit
tests pass, including the1998/1999boundary, missing evidence, deterministic
selection, unknown genes and invalid scores.

Output SHA-256:
`105c3e4b7de08fce7853a787bb25ec50a0e0b5a5b45be7562fa957627da73987`.

```sh
python -m pytest -q tests/unit/test_ob_gpu_batch_witnesses.py
python benchmark_tools/audit_ob_gpu_batch_witnesses.py \
  --root /path/to/orthohmm \
  --output /new/path/ob_gpu_batch_witnesses.json
```

The retained input manifest uses absolute paths; relocation requires the
original inputs at those locations or a separately verified relocation
procedure. The output path must not exist.

## Interpretation Limits

The retained `benchmarks/benchmark_orthobench.py` passes complete species
objects to `search_species_pair` once per directed pair and saves query/
target identities from the returned results. That source is recorded, but
its present contents are not proof of the exact historical executed revision.
Thus exclusion of the defect is conditional on that batching model and
authentic retained outputs. Arbitrary query sub-batches would need their
own witnesses and cannot be cleared from these aggregated rows.

This does not establish historical GPU availability, independently verify
Viterbi scores, distinguish prefilter from scoring rejection, prove absence
of other bugs or clear any other dataset/runtime. It is not a new accuracy
evaluation. No prediction, checkpoint, frozen executor or benchmark score
was modified or rerun.
