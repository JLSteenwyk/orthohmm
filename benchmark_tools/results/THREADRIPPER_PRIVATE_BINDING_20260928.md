# Prospective Private Runtime Binding

The [binding receipt](threadripper_private_binding_20260928.json) combines the
private inference/controller/base snapshot with fresh OS, native-tool and
helper inventories. It covers 57,882 entries and checks 190 distinct frozen
baseline paths, including scientific source, adapter files, tool entrypoints
and OrthoFinder distribution files. The prior private snapshot remains exactly
unchanged. This is a prospective deployment, not historical runtime equality.

Six shared/user roots are explicitly retired: shared Conda bin/lib, terminfo
and OpenSSL configuration, user Python 3.10 site-packages and the user's
matplotlib font cache. Remaining historical inventory roots are retained and
current top-level benchmark Python helpers are included. Retirement is not
proof that no unexercised branch could access those paths.

## Repeated Native Lookup

The [repeated lookup receipt](threadripper_private_lookup_20260928.json) records
two successful executions of the unchanged import inspector using the private
controller Python. Each execution probes both native interpreters under the
new baseline. OrthoHMM loads 913 modules and OrthoFinder loads 1,563. All
observed imported/mapped files and executables are covered by pinned inventory
records with no changed or missing entries. Repeated lookup signatures match
for both interpreters, including module paths, startup hooks, mapping identities
and scientific package origin checks.

Raw inventories, native reports, subprocess output and exact commands remain
local under `benchmarks/work/threadripper_private_binding_20260928`; compact
receipts bind their hashes. The large reports are not a public archive bundle.
Seventeen focused tests pass across root selection, runtime inventory and
native lookup comparison.

## Limits And Next Integration

The native/OS inventory lists 135 symlinks whose targets are outside its roots;
file targets have hashes where available, but directory links are not traversed.
Observed native file coverage is complete only for the declared import probe,
not every workload branch or a complete dynamic-loader closure. This inventory
does not establish a hermetic machine, continuous identity enforcement, quiet
host or controlled timing.

`RuntimeChecker` still launches the inspector with OrthoHMM's interpreter. The
private inference environment intentionally lacks controller-only Biopython,
so integration must use the binding's recorded `controller_python` for that
inspection subprocess. Preserve legacy behavior for old bindings and bind the
changed helper before collector validation. Repeated controller import coverage
and collector v5 native fixtures remain to be checked; do not reuse the old
hard-coded lookup/plan pins in the fixture driver. No production timing or
native scientific inference was launched in this turn; no DGX access occurred.
