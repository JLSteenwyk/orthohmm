# Native TP Losses And Historical Fragment Annotations

[Protocol](NATIVE_QFO_FRAGMENT_PAIR_PROTOCOL_20261006.md) pushed at `b1370417`
before computing the native annotation-stratified counts. Sources and 27 new
tests pushed at `84cd61b1` before actual execution. The
[generated table](native_qfo_fragment_pairs_20261006_v1/TABLE.md),
[full report](native_qfo_fragment_pairs_20261006_v1/report.json),
[10,765-pair ledger](native_qfo_fragment_pairs_20261006_v1/pairs.tsv) and
[independent SQL readback](native_qfo_fragment_pair_readback_20261006_v1.json)
are new descriptive diagnostics, not new scoring/admission or method tuning.

## Observed Changes

Both native P0C0 cells retain initial HMM search on, with profiles/candidate
expansion off. R0 uses group-clique predictions; R1 resolved native pairs.
Annotation bins concern pair endpoints, not all members of a positive family.
Eleven of 563 exact sequence-matched historical records have the original
fragment/incomplete-sequence-positive definition. Both endpoints of every
removed TP pair are unflagged in the historical view.

| Historical Endpoint Bin | Relations | R0 TP | Removed TP | R0 FP | Removed FP |
| --- | ---: | ---: | ---: | ---: | ---: |
| At least one annotation-positive endpoint | 527 | 42 | 0 | 36 | 35 |
| Both matched and unflagged | 10238 | 3081 | 334 | 1672 | 1654 |
| Missing without a positive endpoint | 0 | 0 | 0 | 0 | 0 |

All 42 annotation-positive R0 TPs remain TPs; all 334 TP-to-FN changes occur
among the 3,081 unflagged R0 TPs (10.841%). Of the 1,689 removed FPs, 35 have
a positive endpoint and 1,654 do not. Descriptive FP removal is 97.222% versus
98.923% in those bins. Rates use the appropriate R0 TP or FP denominator,
not every reference pair or an independent-pair sampling model.

This contradicts the specific observed explanation that TP removals principally
involve endpoints with the retained positive annotations. It does not exclude
unannotated biological fragments, establish completeness of unflagged proteins,
validate reconciliation evolutionarily or identify the true cause of losses.
The earlier family-positive bin includes many unflagged pairs: a family-level
positive must not be propagated to every relation.

In the unchanged baseline-only sensitivity view, all 14 later-sequence-version
records become missing. The positive bin stays identical. The other 10,238
relations split into 9,681 unflagged relations (2,960 R0 TP; 323 lost) and
557 missing-without-positive relations (121 R0 TP; 11 lost). All original
relations/decisions remain retained. Zero-denominator rates are null, not zero.
The full report preserves all 96 transition cells and 48 marginal cells,
including zeros. These are not stratified native F1 or causal comparisons.

## Checks And Execution

Primary Counter-based grouping and independent SQLite join/grouping agree
on every ledger row, six view/bin summaries, all transitions/marginals and
12 removal fractions including nulls. Both separately parse the 21,530 native
raw rows and all 563 selected original historical entries, totaling 5,759,359
entry bytes. Each selected entry matches its original hash, accession, taxon,
sequence digest/version, entry version and annotation fields. The two annotation
checks share Bio.SwissProt 1.87; this is not wholly independent parser validation.

Native represented genes, historical membership and original recorded FASTA
identity vectors agree exactly. Prior history/version selection and input
consumption evidence are reused, not reacquired or re-admitted. There are
1,209 original acquisition/source records not newly rechecked by this audit;
their identities remain inherited. Do not describe the selected-entry checks
as a fresh full historical annotation or whole-proteome admission.

The joined [JUnit receipt](native_qfo_fragment_pairs_20261006_v1/native_qfo_fragment_pair_joined_tests_20261006_v1.xml)
reports 119 tests, no failures/errors/skips (3.812s XML; pytest displayed3.84s),
including 27 new cases and existing transition/fragment admission/statistic
contracts. New tests exercise both views, positive-over-missing precedence,
empty bins, null denominators, additions as well as removals, invalid states,
changed truth/universes, canonical-pair duplicates, raw format/member errors
and destination refusal. Earlier 27-case receipt remains retained. Fixture
success is not biological replication; actual runs independently succeed.

Actual review Python3.10.13 uses the retained environment, with Python/site/
preload/library overrides removed and bytecode disabled. The independent
reader uses `-I -B` with Bio.SwissProt available; not stdlib-only or `-S`.
No package installation or scientific-environment modification.
[Primary receipt](native_qfo_fragment_pairs_20261006_v1/native_qfo_fragment_pair_audit_20261006_v1.time.txt)
reports1.17s/62,976KiB maximum process RSS and
[reader receipt](native_qfo_fragment_pairs_20261006_v1/native_qfo_fragment_pair_readback_20261006_v1.time.txt)
1.17s/44,544KiB. Both exit0/zero swaps. These are shared-host postprocessing
costs, not inference costs or an isolated speed comparison. CPU, memory-
bandwidth and I/O contention have unknown, potentially tool-dependent effects.
Fresh capacity was670,187,859,968 available RAM bytes/25,825,280 free swap bytes;
near-full swap and unrelated workloads remain untouched.

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| `audit_native_qfo_fragment_pairs.py` | 13005 | `37afe60dcc6da90a7b4f3c6ef3e7f19116b1e651e5c900a8d8aeca3cadd922b3` |
| `readback_native_qfo_fragment_pairs.py` | 12892 | `85e2c0e98d01ea5fab856959da4fd2c7de7b4e09c8cdd39002dd618a3646f732` |
| `native_qfo_fragment_pairs_20261006_v1/report.json` | 322454 | `36a74ec072917f63d4054f47a84ba0216561a7c1bdbb44c79109a594009f561b` |
| `native_qfo_fragment_pairs_20261006_v1/pairs.tsv` | 739915 | `c38faa7b05e8e575f085bfa31c9486c203aabb0aea0d8968db90e257b5866ac6` |
| `native_qfo_fragment_pair_readback_20261006_v1.json` | 162891 | `d2f6fd3787a1c554ee07f73a981583539fe6a142c9ab19140d70299864614937` |

## Remaining Scope

No new bootstrap draws, pair-IID tests, transferred historical intervals,
independent generalization or biological-completeness claim. These annotations
are development-exposed and may contain errors or temporal limitations. The
result advances native fragment-error and mechanistic evidence (4.3/4.4),
not their full closure, wider valid uncertainty or total HMM contribution.

Failed R1 timing stays ineligible. No original native handle restarted or
unfinished output inspected. Latest check original22444 RUNNING10:45:44;
22445/22450/22451/22452 dependency-pending. Await original terminal gates for
admission/next native identity. Five score cells still unavailable. Original
manuscript/supplement/archive and frozen scientific/execution bytes unchanged;
the full publication goal remains active and completion unproven.
