# OrthoBench Acquisition Support In Source Exports

The optional `orthobench-inputs` source profile closes a concrete packaging
gap: the exported acquisition verifier and rebinder previously lacked their
three required manifest documents. It does not distribute upstream sequences,
reference memberships, dependency binaries or a complete executable study.

## Export And Executed Checks

Source checkpoint `27efba28dfa34898d46bb2411c34b75be68a68ed` adds explicit
schema-v2 selection, three fixed manifest hashes and required helper scripts.
Default source-only selection and historical schema-v1 verification remain.
Build rejects missing/changed support before creating an output tree; offline
verification rejects missing, altered or unexpected payloads and profiles.
The [component guide](../PUBLICATION_SOURCE_COMPONENT.md) supplies copyable
acquisition/rebinding commands using fresh paths outside the immutable bundle.

The three-module focused panel passes **77 tests in 7.38 seconds**, with zero
failures, errors or skips. It exercises both profiles, relocated verification
without Git, default exclusion even when support exists, coordinated manifest
corruption, missing helpers, unchanged revision separation and existing
acquisition/rebinding behavior. The earlier 75-case panel overlaps it.

The [execution receipt](source_orthobench_input_support_20261002.json) records
a real committed export, archive creation, member validation, fresh extraction
outside the checkout and copied-script execution using isolated stdlib Python:

| Check | Result |
| --- | --- |
| Frozen scientific source files | 43, revision `7f3a9e40dd7e79f842cc2c11fb8b548f9a802806` |
| Workflow and support files | 1,778, including exactly three result-document exceptions |
| Indexed payload files | 1,821; 10,552,922 bytes |
| Archive files including index | 1,822; 12,697,600-byte local tar |
| Copied verifier, Git absent from PATH | Successful before and after input preparation |
| Existing upstream checkout | Frozen commit `872d6f30592ab5ff837224db16a514b3f2bb916a` |
| Upstream file verification | 95 files, including 93 benchmark inputs |
| Rebinding | 12 FASTAs, 70 reference files, 11 low-certainty files |
| Ordered basenames, byte identities and nonpath metadata | Unchanged |

Archive SHA-256:
`0d7bc1b0b6da204524f5b7eb84ba92acca9a0e9191057fcc75d9ea8485147628`.
External index SHA-256:
`fa12e1fb49d904208ba646f972f50a2a9c5c8b06e8b37fc17426f83744fc87db`.

The first harness attempt wrote acquisition/rebinding outputs inside the
immutable extracted component. Both preparation stages passed, but the final
inventory guard correctly rejected the extra files. Its logs/archive are
retained. A fresh attempt uses a separate external output directory; no
production guard, scientific input or source byte was loosened or changed.

Separately, the new verifier successfully checks the retained historical
1,707-file schema-v1 component against its original external digest
`f15c343cfb52442d77d1c58ba7292addbe3f190377b7ffe9bdbbd555b122d0b8`.
Two historical invalid-escape SyntaxWarnings remain; syntax compilation is
not import, installation or native execution validation.

## Scope And Remaining Work

This reuses the already acquired upstream checkout. It is not a fresh internet
download, an OS containment test, cross-host restoration, native inference,
new score/uncertainty calculation, public deposition or rights clearance.
The manifests preserve historical plans/results/provenance; their other fields
are not current-table replacements or authorized commands on historical paths.

The [admitted full OrthoBench archive reproduction](RESTORED_ARCHIVE_FULL_OB_RESULT_22377.md)
remains the native evidence and was not rerun. Native wheel/base/installer and
OS dependency closure, full release assembly, unresolved rights and controlled
Threadripper resource measurements remain open. No timing jobs, host-contention
polls, unrelated process/service changes or DGX operations were performed.
