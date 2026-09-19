# Verified Pressure Overhead Deployment

Source/launcher/protocol milestone63d1098 was committed and pushed before
deployment. Exported425 committed files to the fresh remote directory
`/home/jlsteenwyk/projects/orthohmm-publication/pressure_overhead_recipe_v1`.
The transferred archive and local archive match SHA-256
`5d579c315f0ad49abf64dec3e92d328f9b539e0c65d3e0b37a6274b53f0f5ba8`.
Every remote file hash matches its local archive member. No external
symlinks are present. The recipe manifest contains428file/directory records.

- [Recipe manifest](dgx_pressure_overhead_recipe_v1_20260919.json):
  `11ce57a7895d92c456cc2b438157355cb0773cd473f4fbc4fc8305330a33eda3`.
- [Exact-scope authorization](dgx_pressure_overhead_authorization_20260919.json):
  `d1f5878580ed85829b1ce50ba0eebbcfd573975a6eb2b28471762a81978b9dea`.
- [Frozen plan](dgx_pressure_overhead_plan_20260919.json):
  `3950ccfa867c463ccfd7d8693dc331a3c85dc75ba67d9e8db89be4c1dda91c38`.

Authorization covers only the18engineering tasks and explicitly leaves
scientific execution unauthorized. The [remote preflight](dgx_pressure_overhead_preflight_20260919.json)
verified the complete recipe inventory and selected all18tasks using exact
pins; native_execution_started=false. The fresh output root did not exist.
Runtime/environment/input checks remain mandatory before and after each
actual native command, not inferred from this non-executing selection test.

Submit `benchmark_tools/run_dgx_pressure_overhead.sbatch` with `--hold`
and a future eligibility time after preparation. Before release, start a
durable controller-host scheduler recorder from the committed capture
source and verify its live allocation and first retained poll. No DGX
access is permitted during the released panel's quiet window. See the
[prospective protocol](DGX_PRESSURE_OVERHEAD_PROTOCOL_20260919.md).

This receipt does not report inference outcomes or admit scientific timings.
