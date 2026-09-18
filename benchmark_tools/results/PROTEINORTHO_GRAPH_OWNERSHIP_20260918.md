# Proteinortho Graph Ownership Validation

## Result

The original-release native post-clustering graph passes the new complete
input-inventory check:

| Quantity | Validated value |
| --- | ---: |
| Input proteomes | 78 |
| Unique input accessions | 976,504 |
| Species-pair sections | 3,003 |
| Native relation rows | 4,579,157 |

Graph: `qfo_benchmark/results/proteinortho/run/input/qfo.proteinortho-graph`.
Size: 340,068,913 bytes. SHA-256:
`2ecf9df78784fc8add8e99720fb7b07e4f3b5318de7f9d87d0467683e4796cc2`.

All graph and input-file hashes were recorded before traversal and checked
again afterward. Every accession belongs to the species named in its
column's section header. Every relation score is finite. Every unordered
pair of input species has exactly one section; empty sections are allowed.
The existing converter's duplicate and structural checks also pass.

No original converter, graph, prediction or score was changed. A filename
check rejects the usual pre-clustering search-graph path, but filenames
alone do not prove provenance: corrected-run admission must bind this file
to the frozen command, successful execution record and output hashes.
This result does not establish workflow completeness or biological accuracy.

## Reproduction

From the repository root, with the benchmark Python dependencies installed:

```python
from pathlib import Path
from benchmark_tools.validate_proteinortho_graph import validate_graph

print(validate_graph(
    Path("qfo_benchmark/results/proteinortho/run/input/qfo.proteinortho-graph"),
    sorted(Path("qfo_benchmark/input").glob("*.fasta")),
))
```

Focused validation: 17 tests passed across
`test_validate_proteinortho_graph.py` and `test_proteinortho_to_pairwise.py`.
These include species missing from the entire graph, foreign/wrong-species
accessions, nonfinite scores, duplicate sections/relations, malformed rows,
and search-graph rejection. Existing frozen source files remain unchanged.

## Progress Ledger

Completed: Proteinortho complete-inventory validation and historical native
graph check. The preceding FastOMA conversion milestone was pushed as
`bada5d4`.

Scheduler-confirmed running at this turn's check: corrected Proteinortho
`21708`, SonicParanoid `21710`, OrthoHMM `21706_0`; original QfO factorial
assessment `21711` and reconciliation `21671_3`; DGX scaling `21656_14`.
Corrected legacy BLAST `21713` remains pending resources; corrected
OrthoFinder and dependent factorial jobs remain queued.

Next: bind the corrected Proteinortho graph to completed execution and
runtime records, run this inventory check against corrected inputs, then
convert and score. The current historical result does not admit a
corrected-release prediction. Other publication requirements remain open.
