# Full-Native Candidate Trace Replay

## Evidence

The [retained diagnostic](native_candidate_trace_replay_22429/diagnostic.md)
traces job22429/P0C1R0's previously scored58-gene partition difference.
It checks the complete old and native accepted-merge records against the
original63,245-seed partition, round by round. Every endpoint must be an
entire round-start component; conflicting labels and redundant unions fail.
Both complete54,745-group partitions reconstruct exactly, covering all
251,378genes. This is accepted-union reconstruction, not a rerun of search,
candidate eligibility, clustering, scoring or phylogenetic inference.

Both traces have8,500 accepted merges;8,482 semantic merges are common.
Each has six unique round0 merges and12 unique round1 merges. All common
source/target numeric cluster IDs agree. Common accepted hit counts,
coverage and species-overlap values agree;1,725 common support values differ
by at most1.4210854715202004e-14.

Three common round0 anchors have different attachments. Each is at the
four-satellite cap in both traces; the changed-selection support spreads
are0,7.105427357601002e-15 and0. Round1 includes groups with changed
memberships, so not every unique semantic merge implies another final
difference. Only two larger groups and four singletons on each side remain
different at the end. The full58-gene membership difference and every
unique accepted merge are retained in JSON, without selected-example
filtering. No changed gene touches the frozen OrthoBench reference;
reference exclusion and identical aggregate accuracy are inherited from
the [separate score](NATIVE_FACTORIAL_RESULT_22429.md), not rescored here.

## Validation

The JSON has12 direct input/helper pins and SHA256
`c0a5815c10de8ef51e3429af3e08a69bacb89f56ae6386f4fde1f422063f735c`.
A [separate gene-level union readback](native_candidate_trace_replay_22429/independent_readback.json)
uses no replay/partition helper and confirms both whole partitions. It also
checks all12 direct pins and all920 unchanged native helper pins. This is
not a full transitive runtime admission or a new timing observation.

All44 [joined tests](native_candidate_trace_replay_tests_20261004.xml) pass,
zero failures/errors/skips,9.29s. They cover round-start semantics, cycles,
overlap, label conflicts, unsupported score cells, missing/ambiguous trace
bindings, changed checksums, incomplete reconstructions and scored-summary
disagreement, plus actual JSON/Markdown reproduction and the earlier trace
suite. MatchingPython3.10.13 separately reproduces the complete serialized
report exactly. An initial raw Python-object versus loaded-JSON comparison
fails because semantic-key tuples serialize to arrays; comparing at the
JSON boundary passes without rounding or changing scientific values.

Reproduction, from the repository root, with a new output destination:

```bash
env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT \
  PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  benchmarks/work/native_factorial_review_py310_20261004/bin/python -B \
  benchmark_tools/replay_native_candidate_trace.py \
  --preparation benchmark_tools/results/orthobench_factorial_prepared_20260916.json \
  --preparation-sha256 5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382 \
  --score benchmark_tools/results/native_factorial_orthobench_score_22429.json \
  --score-sha256 f1ed516cb1d8f48de35e91aec31384a92d63368e729b5662263018d80e998023 \
  --output benchmarks/work/native_candidate_trace_reproduction_22429
```

## Limits And Next Work

The original seeds suffice to reconstruct both traces; this does not
independently prove the entire fresh seed order. Near ties at capped
attachments are observed, not a causal test of score-bit/order provenance.
Rejected candidates and complete hit arrays remain unreplayed. Do not
change defaults, opportunistically round support or claim a determinism
repair from this evidence. No cause is attributed to CPU contention;
shared-host timing distortion remains unknown and potentially tool-dependent.

Job22430/P0C1R1 remains RUNNING on fresh scheduler evidence. Preserve this
same attempt and all frozen pins. Its actual terminal review and separate
score must precede index4 preparation/release. Other native identities and
the complete uncertainty, generalization and publication package remain
unfinished. This diagnostic neither certifies publication readiness nor
narrows the active goal.
