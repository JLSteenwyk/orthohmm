# Private Controller Integrated With Collector V5

The subsequent [three-path fixture report](THREADRIPPER_PRIVATE_COLLECTOR_PANEL_20260928.md)
records successful satellite_v2 and full OrthoFinder collector checks. The
remaining-scope statements below describe the earlier high-sensitivity milestone.

`RuntimeChecker` now uses a binding's hash-verified `controller_python` when
provided, checking that executable before and after inspection. Historical
bindings without that field keep the original interpreter behavior. The fixture
driver accepts four explicit lookup/plan path-and-hash arguments together and
rejects a plan inconsistent with its runtime binding before creating output.

The [new binding](threadripper_private_binding_v2_20260928.json) and
[lookup receipt](threadripper_private_lookup_v2_20260928.json) preserve the
scientific deployment. Comparing all native/OS/helper entries with the previous
binding finds exactly the [two intended helper changes](threadripper_private_checker_changes_20260928.json).
Private interpreter trees remain unchanged. Two real invocations of the composed
checker pass full tree revalidation and repeated native lookup comparison;
[their receipt](threadripper_private_checker_validation_20260928.json) retains
both independent checks. Twenty-two focused tests pass, including controller
drift before/during inspection, legacy selection and invalid fixture pins.

## Native Collector Fixture

The [submission](threadripper_private_collector_submission_22373.json) launched
one high-sensitivity 16-gene diagnostic, job 22373, through collector v5 with
the private controller and patched inference interpreter. The
[terminal receipt](threadripper_private_collector_22373.json) records COMPLETED
0:0 for allocation, batch and native step. Preparation, before/after runtime
checks and native output validation succeed. Three groups cover all 16 input
genes, with orthogroup bytes identical to the retained historical fixture.

Raw measurement replay passes, including the new native completion boundary,
job-memory observations, reporting-stage evidence, affinity observations and
host-process stream. The completion receipt reports `anchor_only_at_boundaries`.
This is a successful end-to-end collector-v5 integration check on the high-
sensitivity native path, not merely a mocked test or import-only probe.

## Remaining Scope

The fixture ran on a heavily loaded shared host. Its time is not publishable
controlled timing, and no production timing identity was launched. Observed
boundary completion does not prove continuous containment or final job-teardown
memory accounting. The allocation's terminal state is not a substitute for
missing final usage counters. Existing final-accounting and quiet-host gates
remain unmet.

The private-controller v5 fixture has not yet been executed for satellite_v2
or full OrthoFinder. Their earlier native/import checks do not substitute for
that collector integration. Controller optional branches and external symlink
limitations also remain explicit. No unrelated workloads were stopped and no
DGX access occurred. Next, validate the remaining two collector paths and then
address the outstanding production admission requirements.
