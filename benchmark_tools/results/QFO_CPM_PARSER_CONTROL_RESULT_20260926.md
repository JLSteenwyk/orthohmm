# Frozen Parser Controls: Crash Not Reproduced

Both [preplanned parser-only arms](QFO_CPM_PARSER_CONTROL_PROTOCOL_20260926.md)
completed once with exit zero and empty stderr. The [machine-readable result](qfo_cpm_parser_controls_20260926.json)
has SHA-256 `ef73bd0badf9752345d6e7e8a02cb4fa8a1568bac554bd99f30eeefcba589a93`.

| Allocator | Genes / memberships | Groups | Elapsed seconds, descriptive |
| --- | ---: | ---: | ---: |
| Default | 984137 / 984137 | 390845 | 2.0334 |
| Debug | 984137 / 984137 | 390845 | 2.3472 |

Both used the exact frozen `read_partition` function extracted from its
hash-checked module AST, the same retained refined partition and gene universe,
and the retained interpreter/libc identities. Source, protocol, input and runtime
records were checked before/after execution and independently checked again
with both stdout/stderr files. No site initialization, NumPy, BioPython,
igraph, Leiden or OrthoHMM module was imported into either child.

Garbage collection remained enabled at thresholds `[700, 10, 10]`. Both arms
recorded collection-count increases of 508, 46 and 4 for generations 0, 1 and 2,
respectively. Thus successful parsing did not depend on disabling collection.
The counts are observations, not evidence of collection correctness or memory
safety. Each child had a 120-second wall/CPU limit, 4-GiB address-space limit
and disabled core dumps. Ten unit tests cover frozen parser membership checks,
omission of module-level execution, altered inputs, two fixed attempts, timeout
and signal retention, and no-overwrite protection.

```sh
/usr/bin/python3 -S -m benchmark_tools.probe_cpm_partition_parser --root . \
  --output benchmarks/work/qfo_cpm_parser_controls_20260926
python -m pytest -q tests/unit/test_probe_cpm_partition_parser.py
```

Use a new output directory for any separately justified diagnostic; the command
refuses to overwrite this run. This record does not authorize repeated attempts.

## Interpretation and Next Boundary

The exact parser and saved file can complete in these two fresh interpreter
processes. The original crash is not reproduced by this isolated boundary.
This does not reproduce preceding allocation history or scientific imports,
array conversion, native refinement or the original process environment/cwd.
It neither identifies the corrupting operation nor excludes an intermittent
parser/interpreter issue. The AST extraction preserves the function body but
deliberately omits the original module-level imports.

A further diagnostic should isolate the preceding import/native operations
and their interaction with later garbage collection, with intermediate checks
and fixed attempts, rather than repeating full scientific admission until it
passes. The original failed admission 22155 and blocked candidate 22156 remain
unchanged. No optimizer, inference, accuracy scoring, native admission, default
change or DGX access occurred. High-CPM accuracy is still missing.
