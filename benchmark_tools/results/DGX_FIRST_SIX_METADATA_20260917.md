# DGX Completed-Run Metadata Validation

The first six completed scientific runs pass recorded-command and provenance
checks against the frozen 27-run specification. This is not native-output,
resource-replay or scientific timing admission. No runtime rankings are released
from this metadata-only check.

| Array task | Method | Proteomes | Repeat | Metadata | Host evidence |
| --- | --- | ---: | ---: | --- | --- |
| 21656_0 | OrthoHMM high sensitivity | 4 | 0 | Matches | Inconclusive |
| 21656_1 | OrthoHMM phylogeny | 4 | 0 | Matches | Inconclusive |
| 21656_2 | OrthoFinder full | 4 | 0 | Matches | Inconclusive |
| 21656_3 | OrthoHMM phylogeny | 8 | 0 | Matches | Inconclusive |
| 21656_4 | OrthoFinder full | 8 | 0 | Matches | Inconclusive |
| 21656_5 | OrthoHMM high sensitivity | 8 | 0 | Matches | Inconclusive |

Checks include the pinned scientific specification and original command plan;
exact native/GNU-time commands and working directories; original input hashes
and enumeration order; prepared OrthoFinder input copies and post-run copy
identities; pre/post runtime and recipe assertions; successful measurement
status; and requested 20 CPUs, 96 GiB, 85,800-second native timeout, one-second
resource sampling and 30-second host sampling. Preparation and pre/post
verification durations are separate from inference. Live Slurm accounting
confirms completed exit-zero tasks on spark-7ff0 and maps array IDs to the
actual numeric collector job IDs, which are not assumed consecutive.

Only preparation, verification and measurement-summary JSON files were fetched
from these completed tasks. No native outputs or runtime trees were copied or
hashed on the active timing host. These assertions still need raw evidence
replay and native output validation. A successful metadata comparison is not
proof of byte identity for files not independently re-read by this audit.

All six retained host summaries are inconclusive. The existing
[run 00 review](DGX_RUN00_HOST_REVIEW_20260917.md) explains that run's unmatched
kworker-named observations; its explanation is not transferred to the other
five runs without inspecting their evidence. No monitor threshold, classification
or deployed recipe was changed, and no expensive run was restarted or discarded.

Reproduce from the repository root using the downloaded evidence:

```bash
python benchmark_tools/audit_dgx_scientific_metadata.py --root benchmarks/work/dgx_completed_metadata_20260917 --results benchmark_tools/results --indices 0 1 2 3 4 5 --output /tmp/dgx-first-six-metadata.json
python -m pytest tests/unit/test_audit_dgx_scientific_metadata.py -q
```

`dgx_first_six_metadata_20260917.json` records all input hashes, the accounting
snapshot and exact selected run inventory. Fifteen tests exercise each native
method and rejection of changed commands, scheduler/job identities, runtime
hashes, input order/identities, copies, resource settings, timeout state and
invalid durations. Pending or active tasks are neither admitted nor treated as
failures. All 27 runs and the broader publication requirements remain in scope.
