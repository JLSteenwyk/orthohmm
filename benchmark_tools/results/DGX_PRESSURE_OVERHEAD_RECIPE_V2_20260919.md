# Replacement Recipe and Full Preflight

Exported426committed files from17fb1b5 to the fresh remote
`/home/jlsteenwyk/projects/orthohmm-publication/pressure_overhead_recipe_v2`.
Every remote file matches the corresponding local archive member, with no
external symlinks. Archive SHA-256:
`e4754098bb4210d26645e7189f1069413bbf2fff7fd7fdf90b4255f02e479071`.

- [Plan](dgx_pressure_overhead_plan_v2_20260919.json):
  `b644e165dbf4d0beabf1cf4d9b6c314de522e3ebd1b91598ebebea99094c8fff`.
- [Recipe](dgx_pressure_overhead_recipe_v2_20260919.json):
  `50fc3a14c4d5b53c4eb158a22efcfb8b289b346718fad3c78b300172e319d8a1`.
- [Authorization](dgx_pressure_overhead_authorization_v2_20260919.json):
  `17943df08dfe21bf1975ff0dcb29dc66a9507b156e4698b1672e7517cdfd85e6`.

The [non-executing full preflight](dgx_pressure_overhead_preflight_v2_20260919.json)
ran with the exact pinned OrthoHMM-environment interpreter and verified:

- All26673environment/runtime records.
- All10066system-runtime records.
- All429recipe file/directory records.
- Original inputs for each of the three method templates.
- Actual frozen input enumeration and exact four-file order.
- Every task0..17 under the new recipe/plan authorization.
- The output root does not yet exist.

No native inference command was started. The shell launcher explicitly uses
the pinned interpreter; the Python launcher independently rejects a mismatch.
No environment package, native method argument or numerical budget changed.
The old failed recipe/output directory remains intact.

Submit `run_dgx_pressure_overhead_v2.sbatch` held, then verify the durable
controller recorder and first poll before release. Preserve the at-least60s
preparation-to-eligibility delay and remote-access quiet window. Per-run
before/after identity checks remain required despite this preflight.
