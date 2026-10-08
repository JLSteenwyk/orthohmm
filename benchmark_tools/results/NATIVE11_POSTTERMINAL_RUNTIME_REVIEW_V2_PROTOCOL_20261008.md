# Native11 Allocation Aware Runtime Component Successor

The first read-only component invocation, prospectively committed in9cf41015,
exited1 before creating its output directory. The traceback reached
`run_native_factorial_cost.terminal_accounting` and reported
`Accounting identity, resource envelope or terminal state differs`.
Source inspection establishes an import wiring defect: that historical parser
expects `orthohmm_factorial_cost`, whereas native23985 uses
`orthohmm_allocated_factorial`. Both require64 allocation slots; an initial
CPU-envelope interpretation is superseded by this concrete job-name diagnosis.

Retain v1 source SHA256
`5c2e9c80c2c1be1f2972a914c08dcd4669909f755d73b445756f5d3144e6d8e0`
and its tests/protocol unchanged. The separately versioned v2 wrapper imports
the original allocation-aware verifier. It reuses v1 comparison, chronology and
bracket-replay kernels, never v1's defective production wrapper. No source
monkeypatch, gate relaxation, changed scientific setting or native retry occurs.

All original evidence/inventory/lookup and nonadmission requirements in the v1
protocol remain intact. The truthful new schema is
`native11_postterminal_runtime_component_v2`. v2 also creates its separate
destination before preflight so ordinary preflight failures leave failure.json.
The destination is `benchmarks/work/native11_postterminal_runtime_review_20261008_v2`.

Commit and push v2 wrapper/tests/protocol before one new read-only component
execution with the exact source digest. Invoke the existing scientific3.10,
sanitized environment and single native-library thread limits:
`-B -m benchmark_tools.review_native11_postterminal_runtime_v2 --source-sha256 <digest>`.
This is explicit prospective postprocessing recovery of an established wiring
defect, not an automatic inference retry or a new original full-review attempt.
Retain any new failure; do not launch native12 or score from this component alone.
