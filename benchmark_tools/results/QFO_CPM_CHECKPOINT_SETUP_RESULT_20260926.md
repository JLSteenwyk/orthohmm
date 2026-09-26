# Checkpoint/Setup Control: Crash Not Reproduced

The single [preplanned control](QFO_CPM_CHECKPOINT_SETUP_PROTOCOL_20260926.md)
completed with exit zero and all seven boundary markers, ending at
`stopped_before_refinement`. The [saved report](qfo_cpm_checkpoint_setup_20260926.json)
has SHA-256 `794bca9664f0e862ea6bcaac2ff4131af4b54cb77afab3e871c032ffa8109701`.

| Observation | Count |
| --- | ---: |
| Checkpoint genes | 984137 |
| Species | 78 |
| Audited hit rows | 90687327 |
| Saved seed groups indexed | 314274 |
| Saved graph edges mapped | 25501180 |
| Retained refined groups parsed at both probes | 390845 |
| Directed hits selected by broad-panel refinement helper | 0 |

The original frozen numeric auditor verifies checkpoint hashes and numeric
integrity. It reports 983,835 self-hit rows, no nonpositive scores, lexical gene
order and scores between 0.0047840368172403286 and 7.4837209302325585. These
are retained-checkpoint observations, not new sequence-search results or proof
of hit completeness. Both endpoint graph arrays are int32; graph weights are
float64. Species and hit endpoints are int32; hit scores are float64. All seven
arrays remain read-only memory mappings while the final parser runs.

Garbage collection stayed enabled. Generation 0/1/2 collection counts increased
by 1834/166/14 between the first and last markers; no collection was forced.
All 256 checked record entries, child stdout/stderr and five scientific source
records were rechecked after completion. The helper refuses changed plans,
source/control evidence, checkpoint bytes and saved graph/seed inputs. Record
entries include duplicates and are not a claim of a complete import closure.

The elapsed field is 11.7552 seconds for child execution plus parent post-run
identity checks, on the shared host. It is not a controlled inference timing.
The child used the debug allocator, 8-GiB address-space limit, 300-second CPU/
wall limit, original replay cwd/imports and one-thread library settings. No
retry, DGX access, optimization, refinement or accuracy scoring occurred.

Six new tests plus 17 preceding parser/import tests pass (23 total), covering
the exact stop boundary, retained mappings, fixed one-attempt handling of
signals/timeouts/malformed reports, changed plans and no overwrite.

```sh
/usr/bin/python3 -S -m benchmark_tools.probe_cpm_checkpoint_setup --root . \
  --output benchmarks/work/qfo_cpm_checkpoint_setup_20260926
python -m pytest -q tests/unit/test_probe_cpm_checkpoint_setup.py \
  tests/unit/test_probe_cpm_parser_imports.py tests/unit/test_probe_cpm_partition_parser.py
```

## Interpretation and Remaining Boundary

Checkpoint validation/loading, seed indexing, graph mapping and broad-panel
hit selection followed by parsing did not reproduce the failure. This excludes
no intermittent fault and does not reconstruct the original allocation history:
two extra parser probes deliberately occur while numerical objects are alive.
The frozen `refine_cluster_indices` call and writing a newly refined partition
were not executed. Those operations remain the next diagnostic boundary, not
proven causes. In particular, zero selected directed hits does not bypass RBH
graph processing inside refinement.

The original failed admission 22155 and blocked candidate 22156 remain unchanged.
High-CPM accuracy is still missing, with no default change or partial-output
promotion. Success here is not memory-safety evidence or scientific admission.
