# Final Resource Reporting Component

The [actual schema3 archive](final_prepared_resource_reporting_component_20261004.tar.gz)
is built from committed reporting sources at
`a78a5966ead250d20f838d19cce58f851638613e`. It includes the final 27-attempt
table, figure and generated manuscript resource section: 25 measured attempts,
24 eligible observations, six complete cells and exclusions 0/17/20. Three
cells retain missing eligible repeats and unavailable summaries.

Archive: 465,108 bytes, SHA-256
`f48b5a3946bc982f038d65a7a4b2222214c4ef49caacaa8190d48f2a4b11c8a8`.
The [external manifest](final_prepared_resource_reporting_manifest_20261004.json)
is 3,997 bytes, SHA-256
`e36c79c29ed47deec1b32cfded059d92c75726ecc0054fc87c84fc396b551bb8`.
There are 22 payload files plus the manifest, including eight pure reporting
sources. Public archive/manifest match the canonical built bytes exactly.

Restore into a fresh directory outside the repository. Run the copied
`component.py` with the externally retained manifest digest, not a local
collector that follows original evidence paths:

```sh
python -I -S -B /absolute/restored/component.py verify \
  --directory /absolute/restored \
  --manifest-sha256 e36c79c29ed47deec1b32cfded059d92c75726ecc0054fc87c84fc396b551bb8
python -I -B /absolute/restored/component.py replay \
  --directory /absolute/restored \
  --manifest-sha256 e36c79c29ed47deec1b32cfded059d92c75726ecc0054fc87c84fc396b551bb8 \
  --output /absolute/fresh-replay
```

Actual verification uses the reporting interpreter in isolated standard-library
mode from a fresh `/tmp` working directory. Actual replay uses its copied reader
in isolated mode with Matplotlib 3.10.8 and NumPy 2.2.6. It reproduces all three
table files exactly, identical PNG pixels and exact resource-section text.
The [execution receipt](final_prepared_resource_reporting_execution_20261004.json)
records the actual commands, restored reader, replay record and output scope.
Metadata paths are provenance only; replay evaluates copied pure reporting
functions and does not read original evidence paths.

This closes final resource-reporting portability, not raw-accounting or native
inference reproduction. It does not establish isolated speed comparisons,
causal slowdown, final manuscript rendering, a whole-study executable/versioned
release or publication readiness. The shared-host contention caveat and all
scientific limitations remain. Historical interim archives are unchanged.
