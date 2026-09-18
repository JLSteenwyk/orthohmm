# Corrected Proteinortho Runtime and Launcher

Runtime snapshot SHA-256:
`e0fc669cbbd805f9ec768ccd9088977609be3b4a88d9ff95bfb74f5b73f33ae5`.
It binds the frozen command plan, 18 host configuration/helper files,
23 effective container executable/script files and selected effective
container environment variables. Image and runtime executable checksums
remain bound by the command plan. These checks are not a complete host
system-library inventory or historical execution trace.

The runner constructs an explicit environment allowlist, retaining the
prepared PATH/library/thread values and fixing LC_ALL=C. It excludes
uncontrolled inherited container overrides and language startup hooks.
The container reports LANG=C.UTF-8, LC_ALL=C, its own standard PATH and
LD_LIBRARY_PATH=/.singularity.d/libs. This is a prospective reproducibility
control, not proof that every historical environment variable was identical.
Proteinortho's scientific options and native output semantics are unchanged.

Before launch, the runner verifies the command/input provenance and
effective runtime snapshot. It creates a fresh output root and byte-identical
input copies, probes again from the actual native working directory, and
checks the same records after inference. Existing outputs are never resumed
or overwritten. Failure status and logs are preserved. Zero exit status
does not admit native outputs or accuracy.

Validation: 16 planner/runner tests passed, including every runtime-drift
gate, failed native command, explicit environment construction, resource
allocation and overwrite rejection. Shell syntax validation passed. The
real runtime freeze completed successfully. Native inference may be
submitted from the committed pinned executor; conversion and scoring
remain prohibited until independent native-output admission and frozen
conversion/scoring commands. This does not authorize the older destructive
wrapper. The planned run is shared-host, 32 CPUs/192 GiB/72 hours, not
part of the dedicated DGX timing series.
