# Constructor Input-Format Diagnostic

## Evidence And Question

The independently admitted construction panel recorded six native endpoint
mismatches despite intact original and explicit-int64 arrays. Running direct-stage
job21326 has recorded the same mismatches before weight assignment in both minimal
and frozen import modes (partial observations, not a completed/admitted panel).
Neither the optimizer nor profile search was executed in these diagnostics.

The installed Python igraph and bundled C core both report1.0.0; the extension
is _igraph.abi3.so with64-bit integer IDs. Local inspection shows that NumPy
constructor inputs pass through numpy_to_contiguous_memoryview. The upstream
1.0.0 binding has separate memoryview and general iterable conversion branches.
The limited-API memoryview path unfolds to a Python list before conversion;
therefore the proposed alternative is not simply a comparison of direct native
memory versus Python copying. These source observations do not establish which
operation causes the recorded corruption or certify the compiled branch.

Primary source inspected:
[Python constructor](https://github.com/igraph/python-igraph/blob/1.0.0/src/igraph/__init__.py),
[NumPy conversion](https://github.com/igraph/python-igraph/blob/1.0.0/src/igraph/utils.py),
[C edge-list conversion](https://github.com/igraph/python-igraph/blob/1.0.0/src/_igraph/convert.c).
Targeted upstream searches did not establish a matching reported defect. Do not
interpret that search result as proof that none exists, or diagnose a library
bug from a suspicious conditional code path without verifying applicability.

## Fixed Experiment

Six fresh one-CPU workers alternate original NumPy int32 input and a generator
of built-in Python integer pairs, three repeats per format. Both use minimal
imports, the same saved graph, the same order/multiplicity of edges, explicit
vertex count, undirected graph and unchanged weights. Retain loops and isolates.
No sorting, deduplication, resampling, optimizer execution or accuracy selection.

Use the tested direct-stage observer before and after weight assignment, retaining
full mismatch counts, bounded examples, original-array checks and ordered weighted
graph fingerprints. Preserve inputs, sources, library/module identities, actual
affinity, command mode, partial failures and all outputs. The earlier21326 executor
remains frozen and unchanged. No production or historical environment changes.

An iterable format changes allocations and conversion behavior as well as input
representation. Matching workers are not a proven fix, and mismatching workers
do not isolate hardware versus native-library behavior. Before any subsequent
optimizer replay, require full native graph integrity and complete admission;
never select an output because it matches a preferred historical partition.
