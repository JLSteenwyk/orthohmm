# Initial Supplement Assembly Refusal

The first actual assembly invocation at source commit `9888cd1e` exits 1
before creating its output directory. It raises `ValueError: Changed retained
uncertainty protocol` because the retained report's `new_bootstrap_draws` is
an integer count `0`, not the synthetic fixture's boolean `false`.

The 30-test fixture receipt passed but did not exercise that production schema
type. Fix the generator and fixture to require integer zero explicitly; reject
booleans, floats, negative counts and positive draws. Preserve the original
interval report, input selection, scientific results and protocol unchanged.
This is a presentation-layer schema bug, not an interval or scientific failure.

Original invocation uses selection SHA256
`241f9753f9a360fcc0ae84ff14c152f518982104f020d25a0c26f9eaa3237c4a`
and targets `native_qfo_manuscript_supplement_20261006_v1`. No Markdown,
tables or assembly receipt were written. GNU-time failure receipt remains at
`benchmarks/results/native_qfo_supplement_assembly_20261006_v1.time.txt`:
0.23s, 65,624 KiB, exit 1, zero swaps. Initial JUnit remains at
`benchmarks/results/native_qfo_supplement_tests_20261006_v1.xml`.
The original source and fixtures remain recoverable from the committed revision.
Use a new v2 output/measurement path after validation. Do not call the failure
a successful render or repeat any native scientific computation.
