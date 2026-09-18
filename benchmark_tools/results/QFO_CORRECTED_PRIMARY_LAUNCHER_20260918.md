# Corrected QfO Primary Launcher

## Verified Milestone

Prepared fresh-search commands for frozen OrthoHMM high-sensitivity and
OrthoFinder 3.1.5 full inference on the corrected 78-proteome input.
Command manifest SHA-256:
`dbebd3a6915fddeb2b89e0e591a2ac5c798ee6bccec9925fef485260f41baa5a`.
Runtime snapshot SHA-256:
`b619e53430ad5ebd99874961fdfffdf08e3e386b0af5b51cf2bafcdc491c8bd1`.

The runner verifies frozen files, native environment, exact commands and
input inventory before and after execution. It probes OrthoFinder child
tools in the native working directory through the virtualenv entrypoint,
not the resolved system-Python symlink. Current resolution probes are not
historical process execution traces. Output directories must be new;
failed commands and post-run provenance failures remain unadmitted.

Validation: 18 focused tests passed; batch shell syntax passed; real
preflight passed for all 78 input FASTAs, including child-tool resolution.
No corrected inference has been submitted at this milestone. The command
manifest remains explicitly unrun, with execution authorization false.

## Remaining Work

Commit and push the validated launcher, create its pinned executor, and
submit the sequential two-task batch only once. These are shared-host
32-CPU, 192-GiB accuracy runs, not dedicated timing evidence. Independently
validate native outputs before conversion or scoring. All six remaining
comparator rows and eight factorial cells in the corrected-release
protocol remain required; these two runs do not replace that scope.

Original-release jobs 21703 and 21671_2 and dedicated DGX timing task
21656_13 were confirmed running during this continuation. No unrelated
jobs were stopped. TreeFam source retrieval remains unresolved, and its
family-level uncertainty analysis remains unavailable.
