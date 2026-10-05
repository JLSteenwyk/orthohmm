# Full-Native QfO Pair Preparation

## Purpose And Provenance

The older corrected-factorial converters require cached-candidate admissions
and old executor revisions. They cannot be relabeled as admissions of the
forthcoming full-native runs. The new
[adapter](prepare_native_factorial_qfo_pairs.py) accepts only explicit
full-native requests and successful terminal runtime/resource/environment/
semantic-output reviews. It uses a distinct schema and participant namespace:
`ohmm_qfo_full_native_CELL`. No old run or admission is rewritten.

| Native Index | Cell | Submission Meaning |
| --- | --- | --- |
| 6 | P0C0R0 | Cross-species final-cluster clique pairs |
| 7 | P0C0R1 | Native phylogenetically inferred pairs |
| 8 | P0C1R0 | Cross-species final-cluster clique pairs |
| 9 | P0C1R1 | Native phylogenetically inferred pairs |
| 10 | P1C0R1 | Native phylogenetically inferred pairs |
| 11 | P1C1R0 | Cross-species final-cluster clique pairs |
| 12 | P1C1R1 | Native phylogenetically inferred pairs |

P-off retains the initial HMM search. Group clique pairs are not native
resolved ortholog predictions, and root HOGs are not silently expanded into
cliques for an R-on submission. Conversion executes no inference or endpoint
assessment and does not authorize a subsequent inference identity.

## Checks And Retained Results

Validate request/plan/job/index/cell/repeat, all terminal-review categories
and resource scopes, then freshly query successful native scheduler state.
Require the live-controller request comment where available. Check current
helper/output/input hashes and exact agreement with the frozen corrected-QfO
input manifest. Reuse the native output validator's ownership reader to
check original/prepared bytes, complete gene counts, ownership digest and
per-species counts against the independently admitted output.

Require injective UniProt accession normalization for the entire input,
including proteins absent from predicted pairs. Every accession must be in
the frozen reference mapping. R-off reuses the audited group converter and
independent species-size pair count; R-on reuses the strict native pair
writer and independently admitted native row count. Filtered corrected
predictions must preserve every pair and all bytes; mapping losses or count
changes fail explicitly. Both successful and failed conversion receipts are
retained, with no overwrite or automatic retry.

Report submitted pair volume separately from the number/fraction of all
input proteins appearing in any relation. Singleton/no-relation proteins
remain in the coverage denominator. Empty preparations remain explicit;
they are not assigned fabricated zeros for unavailable native endpoints.

Actual preparation requires its own scheduled two-CPU conversion allocation
and a new destination outside the native inference panel. Conversion follows
terminal review and is outside the measured native interval. Its diagnostic
timestamps cover pair writing/filtering/coverage/postflight only, excluding
prior admission and ownership indexing; no standardized whole-conversion
CPU/peak memory or isolated efficiency is claimed.

## Evidence And Remaining Work

[Dependency audit](results/native_qfo_conversion_dependency_audit_20261004.json)
checks the actual seven planned identities against the frozen78-proteome
input hashes, original scoring-environment/mapping hashes and unchanged920
live helper bindings. It does not convert predictions, parse/admit new
full-native outputs, certify fresh scorer runtime or execute native endpoints.
No full-native QfO identity has completed yet.

[735-test report](results/native_qfo_conversion_tests_20261004.xml) passes,
zero failures/errors/skips,11.20s. The55 new tests cover all7 identity
bindings, terminal/resource/semantic refusal, injective/full-input mapping,
group versus native predictions, independent counts and protein coverage,
empty predictions, preservation of partial failures, joined synthetic
conversions, scheduler/allocation refusal, mapping losses and no-overwrite.
Joined scheduler/plan tests use explicit doubles, not actual native admission.
One initial empty-fixture filename collided with its own source; retain that
1failure/44pass report in work and use a distinct test destination.

After actual QfO inference/review/conversion, add a separately bound native
six-endpoint assessment and independent score admission for this new schema;
the old cached-assessment gates must not be bypassed or relabeled. Check all
scorer/reference/runtime dependencies and existing FAS population conventions
at that handoff. GO/EC similarity and FAS are not F1; the six-endpoint mean is
project-defined and secondary. Native scoring, appropriate uncertainty and
the full publication requirements remain uncompleted by this adapter.
