# Staged Native OrthoMCL Inputs

OrthoMCL 1.4 mode4 derives the index names from the BPO filename and uses
existing indexes without rebuilding them. The parallel patch stores pair
caches beside that BPO. Therefore production inference needs a fresh,
isolated directory with exactly the admitted BPO, its offset/query indexes
and the original admitted GG mapping, not a reused search/cache directory.

`stage_orthomcl_native_inputs.py` copies these inputs to the native layout:
`all.bpo`, `all_bpo.idx`, `all_bpo.se`, `all.gg`. It validates source records,
checks each copy before publishing its filename, checks gene/species scope,
and reruns the complete native index validator on the copied files. Original
and staged hashes are checked again afterwards. It rejects existing output
directories, dangling destination symlinks, shell-unsafe paths, invalid
scope, missing/aliased sources, copy corruption and index disagreement.
Failures and partial files are retained, with no implicit resume.

Files are separate copies, not hardlinks or symlinks. The local filesystem
had about12TiB available when checked; actual production capacity still
needs checking against completed artifact sizes before staging.

## Actual Native Fixture

`probe_orthomcl_staged_inference.py` binds the source fixture to the retained
three-arm report, verifies dedicated Python and native runtime inventories,
stages its prebuilt indexes, and runs guarded native mode4 with two pair
workers in a fresh tool directory.

The installed-runtime probe passed:

-42input proteins across3species;379BPO records,42query ranges and380offset
 entries including EOF;24,168BPO bytes.
-12final groups containing41proteins, exactly matching the untouched serial
 fixture partition. The remaining input protein is not silently added.
-Every staged input/index retained identical bytes and modification time
 after native inference, supporting that the supplied indexes were not
 rebuilt during this run.
-Python and Perl/system runtime checks passed before and after. The native
 execution exited zero. The earlier bundled-example partition discrepancy
 remains documented, not resolved or hidden by this fixture.

61focused tests passed with legacy-runtime tests enabled, including a real
dedicated-interpreter staging/inference test and staging failure paths.
These tests do not prove all64-worker schedules or full-data accuracy.

Retained report: `orthomcl_staged_native_inference_probe_20260918.json`,
SHA-256 `2f83d275bedae2052b7341035a7c4a2c4a18cd174e6a6aaf11ae489de53b262d`.
Final raw artifacts are in
`benchmarks/work/orthomcl_staged_native_inference_probe_v2_20260918/`;
the earlier probe directory is preserved.

## Remaining Production Work

The staging component requires admitted source records and runtime checks
from its caller; it does not infer that arbitrary input records are admitted.
Bind it to successful checkpoint admission21749, the frozen configured
native source manifest and scheduler resources in the full inference
wrapper. Then independently validate final native groups and run conversion
and QfO scoring. No production inference was launched by this fixture work.
