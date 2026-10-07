# Complete Candidate Pair Paths With Original-Protein Alias Proof

The [prespecified follow-up](NATIVE_QFO_CANDIDATE_ALIAS_GROUP_PROTOCOL_20261006.md)
resolves the initial diagnostic's exact identifier boundary using original
QfO protein identity, not suffix guessing. Protocol9102ed47 is pushed before
selected alias outcomes; tested exporter/independent reader173a4091 is pushed
before the actual new analysis.49new tests and83prior grouping/diagnosis tests
pass132/3.14s before execution. The original failed exportbfad1eef and prior
diagnosis3bc6646e remain unchanged; neither namespace nor source is retried or
overwritten. This is an explicitly versioned, evidence-backed diagnostic join.

## Identity Proof

Both original conversion inventories bind the same2020 mapping1c10f6ce.
The candidate databasecd938ce7 is the original admission/execution-inventoried
artifact. Original map_relations/VGNC sources remain originally inventoried.
The two exact scored aliases and uniquely matching normalized native genes
share protein number and species in BOTH original map and retained SQLite
rows; the complete native universe has exactly one accession per selected
protein number. No name-derived substitutions, current online annotations or
queries of predicted ortholog edges are used.

| Scored Accession | Native Accession | Original QfO Protein | Species | Original DB Row IDs |
| --- | --- | ---: | --- | --- |
| Q17QN5_BOVIN | Q17QN5 | 594,577 | BOVIN | 698207, 698208 |
| Q1RMT5_BOVIN | Q1RMT5 | 594,851 | BOVIN | 698557, 698558 |

Native genes are`tr|Q17QN5|Q17QN5_BOVIN` and`tr|Q1RMT5|Q1RMT5_BOVIN`.
All selected SQL rows, exact bridges, original source/input records and false
rescore/admission/readiness flags are in the report. This proves equivalence
within the retained QfO identifier system, not independent biological validity
of that reference mapping or a general alias-resolution method.

## All Changed Pair Paths

| First Connected Round | Path | Recovered TP Pairs | Added FP Pairs |
| ---: | --- | ---: | ---: |
| 0 | Direct cross-endpoint | 131 | 1,794 |
| 0 | Transitive union | 11 | 98 |
| 1 | Direct cross-endpoint | 19 | 239 |
| 1 | Transitive union | 1 | 2 |
| Total | All paths | 162 | 2,133 |

The complete353638-group partition reconstructs from394328baseline groups
and40690accepted unions with all984137genes retained. The union table's all
42080rows and ALL2295changed pairs are checked, with no omitted/mapped subset,
self-collapse or alias-collapse. Original scored pair IDs/categories remain
unchanged in the full ledger alongside native genes, original baseline/group
keys, first connected round, direct event index and path. Round indices0/1
are software rounds, not bootstrap replicates or phylogenetic times.

Of2,133added scored FPs,2,033are direct and100transitive; of162recovered
TPs,150are direct and12transitive. Therefore an explanation in which all
added scored FPs occur only through transitive attachment is inconsistent
with these paths. A direct cross-endpoint GROUP relation is not necessarily
a direct sequence/HMM hit for that particular protein pair. Transitive paths
do not establish absence of search evidence. Candidate eligibility, rejected
alternatives and numeric support are not recomputed; no support mechanism,
true homology/duplication, calibrated confidence or causal accuracy gain is
claimed. The original precision/recall/F1 trade-off remains unchanged.

## Independent Readback

The separate stdlib reader derives bridges again from original mapping/SQL
rows, parses groups separately and traverses an undirected group graph rather
than importing the exporter, original primary or accepted-union kernel.
It verifies all2295ledger rows and all eight state/round/path summary cells,
whole partitions, original inventories and55input/source records before and
after reading. Original scored orientations are preserved even if canonical
native accession order reverses. Its actual execution exits0 and agrees
exactly, not merely on pooled counts. No raw scorer, conversion or transitive
accuracy admission is repeated.

## Identities

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| New alias-join report | 24,497 | `3b3ab8f1132e70362e1450ac990caf7a9dad536e528bc54d6af39db16183cd53` |
| Complete changed-pair ledger | 391,881 | `3f0598f7e00cfc0942146604ab1863ad9f03f2d5017efd5244c7e76a29d43642` |
| Independent readback | 3,942 | `4afc4aa6853a69027b9dc710602309f44923315b1e68dbae4cd1e61866f4b92d` |

[Report](native_qfo_candidate_alias_group_20261006_v1/report.json),
[complete ledger](native_qfo_candidate_alias_group_20261006_v1/changed_pair_groups.tsv)
and [independent reader output](native_qfo_candidate_alias_group_readback_20261006_v1.json)
are immutable. Primary32.26s/2632040KiB and reader22.35s/2631756KiB each
exit0/zero swaps. These are shared-host postprocessing observations, not
tool-speed evidence, native resource costs or isolated efficiency rankings.
Available memory659396136KiB and full swap were observed before launch;
the existing inference and unrelated jobs are untouched.

## Reproduction

In fresh paths only, use the retained scientific Python environment; do not
rerun or overwrite the original failed export. After the new primary finishes,
pass its actual report SHA256 to the separate reader. No package install is
required. For either command, use the common sanitized environment:

```bash
env -u PYTHONPATH -u PYTHONHOME -u PYTHONUSERBASE \
  -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -I -B \
  benchmark_tools/join_native_qfo_candidate_aliases.py \
  --diagnosis benchmark_tools/results/native_qfo_candidate_group_id_diagnosis_20261006_v1/report.json \
  --diagnosis-sha256 3bc6646e1dc7cb6488ceb6aaaf763e3a4e3cdcbd41268f2066ef2d8008f5be74 \
  --output NEW_FRESH_ALIAS_JOIN_DIRECTORY
```

Use`-I -S -B`for the stdlib reader with`--report`,`--report-sha256`and a fresh
`--output`file. Tests use the existing test environment, not the scientific
environment without pytest. Neither replication command submits inference,
changes metrics, reruns scoring, fixes failed timing or establishes valid
VGNC confidence intervals/independent generalization/publication readiness.
The full seven-part goal remains incomplete and active.

Final178joined tests pass12.14s with zero failures/errors/skips, including
actual source/input/ledger/readback identities and exact manuscript/result
table correspondence. Bound analysis-source bytes remain unchanged after
execution. Prior failure/diagnosis contracts now describe their INITIAL
limitation separately from this completed, evidence-backed follow-up.
