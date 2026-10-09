# Prospective Fragment Runtime Inventory Amendment

The first current-runtime preflight against the retained September manifest
exited 1 before inference with `Package inventory changed: orthohmm`.
The frozen source, native-library, profile-probe and executable-file checks
preceding that inventory comparison passed. Preserve this refusal; no
fragment inference ran, and this is not an inference retry or a changed
historical admission.

A subsequent bounded inventory diagnosis found the same Python versions for
both tools and an identical full OrthoFinder package inventory. The OrthoHMM
inventory has added, removed and changed distributions relative to September.
It must not be labeled historically identical. Do not install or remove
packages in the shared environment to recreate obsolete notebook/tooling
dependencies, and do not rewrite the original manifest.

The separately named prospective `controlled_fragment_execution_runtime_v1`
manifest may adopt the measured current OrthoHMM inventory only if Python and
all explicitly required scientific distributions are unchanged: NumPy,
Numba, llvmlite, Biopython, SciPy, igraph, leidenalg, DendroPy, psutil,
texttable and packaging. Record every inventory difference, including absent
distributions. OrthoFinder's whole inventory must still match exactly.
These checks are enforced by the tested controlled_fragment_runtime.py.

All frozen core source, native binary, external executable and distribution
identities, commands and scientific settings remain unchanged. The existing
strict environment verifier must pass against the explicitly new inventory,
including the native profile probe and executable resolution, before any
launch. The new runtime is for the prospective fragment condition only; it
does not retroactively reinterpret baseline environment equality or establish
historical output equivalence. The paired analysis must disclose that limit.

This amendment addresses an observed input/runtime contract change rather
than adding an all-tool hermeticity or operating-system gate. No host service,
unrelated workload, frozen scientific source or completed result is modified.
