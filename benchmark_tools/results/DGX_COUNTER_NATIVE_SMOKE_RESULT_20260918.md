# Counter-Based Native Smoke Results

Array21802 executed three sequential exclusive20CPU/96GiB tasks on spark-7ff0.
All completed0:0 with zero restarts. Code/protocol f2904ea and batch entrypoint
bdd9735 were committed and pushed before submission. All14transferred recipe
files were compared with committed bytes before execution. Recipe manifest
SHA256:9b4475bf95f0cdbf54f2c2cb56238d194c5db59a0da92a9ef7b63739c63ecaab.

| Configuration | Input proteins | Groups | Root groups | Native pair rows | Observer snapshots | Counter errors |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| OrthoHMM high-sensitivity | 645 | 98 | Not applicable | Not evaluated | 5 | 0 |
| OrthoHMM satellite_v2, inferred phylogeny | 645 | 98 | 98 | 1835 | 8 | 0 |
| OrthoFinder3.1.5 full | 645 | 99 checkpoint groups | Not evaluated | 1834 | 11 | 0 |

These are validity/count observations, not accuracy scores. The native command
wall observations3.090376294,6.925679675,9.790886231seconds are engineering
diagnostics only. They are not a speed ranking, scaling measurement or proof
of negligible observer overhead. No timing ratios or significance tests were
computed. Native step and observer batch cgroups were distinct; both shared
the allocation's CPUs. All before/after runtime and original-input identities
and native input enumeration matched. OrthoFinder's fresh input copies were
independently checked against original basenames and bytes.

## Validation

`audit_counter_native_smokes.py` checked the exact relocated commands against
the pinned original smoke specification, terminal scheduler status, recorded
runtime checks, raw host summary replay, scope separation, work bracketing,
and the existing native group/pair/graph and GNU-time validators. Both
OrthoHMM partitions cover the exact645-protein universe; satellite root groups
also cover it, and native pairs match the native metrics. OrthoFinder3.1.5
reported completion and passed checkpoint, input-copy and finite graph checks.
Tree equivalence and prediction accuracy were not newly established.

The first audit attempt incorrectly required identical preflight/postflight
keys for OrthoFinder. Its fresh copied-input evidence exists only after
preparation. The audit was corrected to compare common runtime/original-input
identity fields and separately validate copies through the native validator.
No job was rerun, no raw result changed, and the failed attempt wrote no
success report.

Machine-readable evidence is `dgx_counter_native_smokes_21802.json`; it retains
the scheduler text, verification/counter reports and native validation records.
`dgx_counter_native_recipe_20260918.json` records the remote recipe inventory.
Raw archive: `benchmarks/work/dgx_counter_native_smoke_21802/`,1063files,
3820668bytes at audit. Raw observer streams and complete native outputs remain
there. Seventy-eight focused tests pass, including retained result replay.

## Limits And Next Step

This completes a functional counter-monitor smoke on all three target
pipelines with unchanged scientific flags. It does not establish quiet-host
admission, counter-window error bounds, general overhead, exact RSS, or
comparability of long runs. The original27timings remain descriptive.
Before any new scaling panel, freeze a prospective timing inclusion and
execution protocol that addresses these remaining limitations. No such
panel is authorized or submitted by this smoke. Publication readiness and
controlled-timing admission remain false.
