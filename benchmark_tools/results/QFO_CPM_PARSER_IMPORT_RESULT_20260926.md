# Import-Context Parser Controls: Crash Not Reproduced

Both arms in the [fixed protocol](QFO_CPM_PARSER_IMPORT_PROTOCOL_20260926.md)
completed once with exit zero. The [saved result](qfo_cpm_parser_imports_20260926.json)
has SHA-256 `80ca73e6a9557a00c51940c3778ae0765bf9f8db19891f09175de34f5f2788b9`.

| Context | Genes / memberships | Groups | Elapsed seconds, descriptive |
| --- | ---: | ---: | ---: |
| Site initialization, no explicit scientific imports | 984137 / 984137 | 390845 | 2.2439 |
| Site initialization and frozen replay imports | 984137 / 984137 | 390845 | 2.7120 |

Both use the debug allocator, original replay cwd/PYTHONPATH, one-thread
numerical settings and the same frozen parser helper, universe and retained
partition as the earlier no-site controls. Each has a 120-second wall/CPU
limit and 4-GiB address-space limit. All three markers (`before_imports`,
`before_parser`, `after_parser`) appear in both stderr files, with no other
stderr content. No retry occurred.

The site-only inventory includes `site` but not NumPy, BioPython, OrthoHMM,
igraph or Leiden. The frozen-import arm includes `site`, NumPy, BioPython and
OrthoHMM, with five explicitly resolved scientific source records under the
frozen replay checkout. Neither arm imports igraph or Leiden. This is the
refinement worker's import sequence, not the separate clustering worker's.

Garbage collection stayed enabled at thresholds `[700, 10, 10]`. Generation
0/1/2 collection counts increased by 508/46/4 in each arm. All 237 checked
record entries were revalidated, along with child stdout/stderr and the five
scientific source records. Duplicate identities may occur in the evidence
list; 237 is not a count of distinct files or full transitive runtime coverage.
Post-run comparison also confirmed that all five resolved scientific source
records exactly match their entries in the original Memcheck diagnostic's
checked records, and the frozen checkout has no tracked source differences
from its HEAD. The runner's preflight checks core launcher files and recorded
site-package files; this additional five-module historical comparison was a
separate post-run check, not a claim of exhaustive import-closure preflight.

Seventeen tests pass across these controls and the preceding parser-only
controls, including fixed arms, limits, changed evidence, timeout/signal
retention, malformed completion output, partition validation and no overwrite.

```sh
/usr/bin/python3 -S -m benchmark_tools.probe_cpm_parser_imports --root . \
  --output benchmarks/work/qfo_cpm_parser_imports_20260926
python -m pytest -q tests/unit/test_probe_cpm_parser_imports.py \
  tests/unit/test_probe_cpm_partition_parser.py
```

## Bounded Interpretation

Importing the frozen refinement context and then parsing the saved partition
did not reproduce the crash in this observation. This does not exclude an
intermittent fault, reproduce checkpoint allocation history, establish memory
safety or admit high-CPM scores. Numerical checkpoint loading, array conversion
and graph-based refinement remain untested by these two arms.

Source inspection of frozen `refine_cluster_indices` shows NumPy/Python graph
processing, without direct igraph/Leiden calls. Its broad-dataset branch still
processes RBH graph arrays when `production_refinement_hits` returns empty
lists; zero directed refinement hits must not be interpreted as no refinement
work. A further bounded diagnostic should isolate checkpoint/array operations
and intermediate refinement boundaries with retained observations, rather than
repeat full admission until it succeeds. Original admission 22155 remains failed
and candidate 22156 remains blocked. No DGX access or scientific changes occurred.
