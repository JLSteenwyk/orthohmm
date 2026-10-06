# Native SwissTrees Observed Graph Support

## Actual Finding

The complete 2,023 changed P0/C0/R0-versus-R1 SwissTrees pairs were traced
through the retained final graph files, without changing the method or
rerunning inference. The [report](native_qfo_swiss_graph_support_20261006_v1/report.json)
and [complete case TSV](native_qfo_swiss_graph_support_20261006_v1/pairs.tsv)
preserve both native views, every selected graph record and shortest-path
witness, all 23 original candidate families and their 1,139 native members.
These are native ablation cells with initial HMM search on, not relabeled
selected-default/high-sensitivity competitor results.

Each complete graph has 25,630,303 canonical undirected edges, matching its
originally admitted final network-edge metric. Both retained files have
identical current size and SHA256, and each contains 23,235 induced edges
within the selected candidate families. All 51,260,606 original graph rows
were scanned by both implementations. No nonpositive weights were found;
this was an observation, not an exclusion criterion. Original file bytes
were not rewritten or replaced.

The observed cross-tabulation is identical in the two native views:

| Removed Prediction | Direct Significant Search Support | Direct Final Graph Edge | Indirect Within-Family Path | No Within-Family Path |
|---|---|---:|---:|---:|
| TP | None | 0 | 60 | 0 |
| TP | One direction | 6 | 0 | 0 |
| TP | Both directions | 264 | 4 | 0 |
| FP | None | 0 | 479 | 0 |
| FP | One direction | 16 | 4 | 0 |
| FP | Both directions | 1,121 | 69 | 0 |

Thus 270/334 removed TP and 1,137/1,689 removed FP have a direct final graph
edge. The remaining 64 TP and 552 FP have indirect paths within their actual
original candidate families. All 539 pairs lacking a significant direct hit
have such paths; absence of a direct hit alone does not establish an
unsupported assignment or isolate a prefilter/scoring failure. A further
77 pairs have significant direct hits but no final direct graph edge.
The final combined graph does not isolate which reciprocal-best-normalized-hit
or singleton-attachment rule accounts for those individual edge absences.

Of all changed pairs, 1,407 have shortest distance one, 607 have distance two
and nine have distance three in each view. These are unweighted induced
graph distances, not evolutionary distances, orthology confidence or accuracy
scores. Complete native family members are retained as potential intermediates,
including same-species and non-reference genes; the traversal does not restrict
bridges to the changed/reference endpoints. No true duplication history or
biological correctness is inferred from this connectivity.

## Evidence Scope

The final graph files were not inventoried by either original output validator.
They are newly observed current-byte evidence, not retrospectively admitted
original graph files or proof of continuous integrity during execution.
They are checked against original metrics, a valid bound gene universe and
the original candidate partition. Their current agreement strengthens stage
localization, but must not be promoted into a claim that graph bytes were
known unchanged throughout the original runs.

The admitted R0 partition has 394,328 families and 984,137 unique genes. Its
original SHA256 matches the already reconstructed R1 phylogeny-input partition.
Every selected case's exact native gene names, family identity and original
TP-to-FN or FP-to-TN transition match the bound search/localization reports.
All graph weights on direct changed-pair edges match an actual stored
significant-search score. Those scores remain homology support, not calibrated
orthology confidence. This diagnostic does not re-read the 181-million-row
significant-hit arrays or repeat their previously verified extraction.

Graph presence and connectivity are not causal tests of cluster recruitment,
Leiden/refinement behavior, tree correctness or reconciliation biology. The
final graph precedes refinement and combines normalized-hit and singleton
edges; it is not an isolated initial normalized-hit graph. The analysis is
retrospective and development-exposed, with dependent pairs. It adds no
accuracy endpoint, uncertainty draw, subgroup significance or generalization
claim. A missing induced path would not establish global disconnection;
none occurred in this cohort.

## Verification

The [post-hoc protocol](NATIVE_QFO_SWISS_GRAPH_SUPPORT_PROTOCOL_20261006.md)
was committed as `fbeb4f4a` before inspecting selected graph outcomes.
The implementation and 49 passing synthetic/contract tests were committed
and pushed as `6d9643d2` before either actual scan. These fixtures are not
native scientific or resource admission.

The primary checks complete canonical row order, duplicates/self edges,
known endpoints, finite weights and exact original metric counts. It extracts
all selected within-family edges, performs deterministic breadth-first search
and saves full path witnesses, row offsets and both native views. Relevant
private pipeline/graph/helper sources match the frozen plan and deployment
baseline. This is a direct, partial source/evidence check, not re-admission of
all frozen helpers, input proteomes or original execution.

The [independent readback](native_qfo_swiss_graph_floyd_readback_20261006.json)
imports neither primary scanner nor traversal. A separate CSV scan checks
every complete graph row and all 46,470 selected induced-edge records across
both views. NumPy Floyd-Warshall matrices independently establish every case's
shortest distance; each submitted path must be simple, have the correct endpoints,
use actual edges and attain that distance. Matrix allocation is bounded to
2,048 members; the largest actual selected family has 189. The reader checks
all 2,023 cases, 4,046 TSV rows/headers, 36 cross-tab cells including zeros,
source bindings, member inventories and graph-agreement scope. It verifies
direct evidence digests before and after its scan.

The final joined suite passes all 444 tests in 20.80s, without failures, errors or
skips. New graph tests contribute 50 cases, covering malformed full graphs,
duplicates/unknown endpoints, nonpositive weights, candidate/graph/source
bindings, direct-hit weight requirements, deterministic paths, non-reference
same-species fixture bridges, allocation bounds, a complete source-bound
synthetic 2,023-case export/readback and deliberate corruption. The actual
result/receipt test checks retained source/report identities and counts
without repeating another large graph scan. Existing mechanism figure,
sequence strata, search, reconciliation, transition, count, uncertainty-binding
and scientific export contracts remain joined. Retain the earlier 40- and
49-case XMLs and the earlier 444-case/22.23s receipt; no failed actual graph
attempt occurred. The final suite replaces an untracked-private-directory
dependency in synthetic fixtures with generated inert source-identity files.
Those files are never executed. This improves fixture reproducibility without
changing scanner sources, actual results or repeating scientific scans.

Actual primary and independent scans use sanitized Python 3.10.13; independent
NumPy is 2.2.6. GNU-time receipts record 86.09s/282,648 KiB and
112.16s/355,140 KiB, respectively, with zero swaps and exit 0. These are
shared-host postprocessing observations, not native inference timing or an
algorithm/tool speed comparison. CPU, memory-bandwidth and I/O contention have
unknown, potentially tool-dependent effects. Roughly 627 GiB RAM and 9.5 TiB
disk were available before launch; nearly full host swap remains disclosed.
Failed R1 inference timing stays ineligible. No unrelated workload changed.

Primary source SHA256:
`ec2370581636796645bf8fdd9300e093d689c4fc1ab97682ad8cf0198c5b9894`.
Independent source SHA256:
`1b8176e2e3a79d6d6e3656cfe9eca88498d02e315d904056a5fb260de8d8a354`.
Each observed graph: 1,724,697,831 bytes, SHA256
`4f5d3d0ab83e9bb39a8e53c18d9af2ddcbf9ef508caf73631333b2a343c28e6e`.
Report: 12,858,109 bytes, SHA256
`014455a42df4c60e70025d763e1d7920f842bb55a67b8ac6697100cdda7e9ca0`.
Case TSV: 602,392 bytes, SHA256
`aa4fc300907201d674b1bce1bff8361861322571f5986ae3ad11147298c90c13`.
Independent readback: 14,484 bytes, SHA256
`5fba0de77a8c609dcb8b872f7972b77c0e4b95e2c5a0a43ef51f68bdc99f760d`.
The report includes selected diagnostic data, not the full 1.7-GB graph files.
Full original raw graphs are retained locally and are not committed.

## Reproduction And Remaining Work

From the repository root, with the original bound local artifacts:

```bash
env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  benchmark_tools/trace_native_qfo_swiss_graph_support.py \
  --search benchmark_tools/results/native_qfo_swiss_search_support_20261006_v1.json \
  --search-sha256 787bde01ffacd60547fb290bdbb1d04242230ed05ce77d25d270227310252d5a \
  --readback benchmark_tools/results/native_qfo_swiss_search_support_code_readback_20261006.json \
  --readback-sha256 8ffc58fe192fa060de4e068161e72b8ab626066446608e5c0011925e62c67ede \
  --plan benchmark_tools/results/native_factorial_receipt_amendment_20261004/plan.json \
  --plan-sha256 6c87babcbb5581830e0b9e7b9bf9aaba30a85bde1c4ab465e561017e67e9c89b \
  --output benchmark_tools/results/native_qfo_swiss_graph_support_replay

env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  benchmark_tools/readback_native_qfo_swiss_graph_support.py \
  --report benchmark_tools/results/native_qfo_swiss_graph_support_20261006_v1/report.json \
  --report-sha256 014455a42df4c60e70025d763e1d7920f842bb55a67b8ac6697100cdda7e9ca0 \
  --output benchmark_tools/results/native_qfo_swiss_graph_readback_replay.json
```

Use fresh paths. The second command checks the retained result, not the new
replay; supply the new report and its actual digest to check that instead.
Local absolute evidence paths/private source copies remain requirements;
this is not portable whole-study archive restoration. Original manuscript,
PDF and archived components are unchanged and do not already include this
companion evidence. Update complete assembly and archive only when warranted,
not by relabeling an earlier component as containing new evidence.

Original 22444 was RUNNING at 7:44:12 and 22445/22450/22451/22452 remained
dependency-pending. No unfinished output was read, job restarted, successor
released, scientific source/default/endpoint changed or new native admission
made. Remaining native cells, matched-search interactions, tree/other error
strata, valid wider uncertainty, independent generalization, original TreeFam
files, provenance, complete manuscript/reproducibility/release/deposition
requirements remain open. Shared-host contention is authorized, not a
quiet-window or dedicated-host blocker. The full publication goal remains
active with completion unproven.
