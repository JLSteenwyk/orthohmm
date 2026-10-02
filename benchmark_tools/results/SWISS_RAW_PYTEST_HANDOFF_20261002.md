# Raw SwissTrees Pytest Handoff

The [raw relocation workflow](SWISS_RAW_SOURCE_RELOCATION_20261002.md) already
exports tables from exact restored files. Extend its handoff to the four
retained raw-export regression cases, preserving their existing statistics,
row-equality, wrong-count-hash and no-overwrite assertions. No production,
scientific, raw manifest, counts, scores or historical archive bytes change.

## Explicit Options

Add optional path/digest pairs for each raw panel. A non-autouse fixture forwards
only explicitly supplied arguments to the duplication/fragment exporters.
Descriptive/identity tests receive no extra arguments. Without options, raw
tests keep their original path checks; no automatic lookup, hash recomputation,
source substitution, skip or mock is introduced. Partial pairs are usage errors
before collection. Production still validates the digest-pinned binding, complete
ordered original records, actual restored files and source/statistic identities.
Configuration forwarding alone does not validate the files. Binding integrity
uses an independently retained SHA256, not a digital signature; `signed` below
is only the recorded test configuration's label.

Use the complete transferred directories from the relocation workflow and
their independently retained binding digests, not digests calculated from an
untrusted received manifest. From a restored checkout:

```bash
python -m pytest -q \
  tests/test_export_swiss_duplication_strata.py \
  tests/unit/test_swiss_strata_count_selection.py \
  --swiss-duplication-bindings /restored/duplication-inputs/bindings.json \
  --swiss-duplication-bindings-sha256 2009be32157feddd2d0f6f4191bcbf5e0d33416b24e9a50e402f6aa023df8e58 \
  --swiss-fragment-bindings /restored/fragment-inputs/bindings.json \
  --swiss-fragment-bindings-sha256 0319644fbb580d87ab123676a7220f89f053cac89e75a5674c4295dd4d3f72f1
```

These digests identify the exact private directories retained in
`benchmarks/work/swiss_raw_sources_{duplication,fragment}_20261002`, not publicly
available downloads. Do not silently use them for a newly generated, different
binding. No raw data/native image is committed, uploaded or granted new rights.

## Verification

Thirteen configuration tests cover defaults, unchanged forwarding, independent
panels, preserved markers, partial pairs and unknown exporter kinds. **27 default
local tests pass in 3.69s**. The configured nine-module panel has **200 passes
in 11.42s**, zero failures/errors/skips, and overlaps that earlier selection.
Existing synthetic fragment tests intentionally do not consume the raw fixture.

Stage 1,825 source/data files (12,937,857 bytes) and copy both existing raw
directories into a temporary tree. Four fresh isolated Python 3.10.13 children
test distinct configurations; retain all failures and process statuses:

| Configuration | Actual pytest outcome | Forbidden opens after canaries |
| --- | --- | ---: |
| No bindings | 23 pass/four fail in 0.65s; exit 1 | Four blocked original-path attempts |
| Correct digest-pinned bindings (`signed`) | 27 pass in 7.16s; exit 0 | Zero |
| Wrong binding digest | Selected export fails in 0.12s; exit 1 | Zero |
| Incomplete option pair | Usage error before collection; zero tests; exit 4 | Zero |

The unconfigured failures are exactly the four retained raw cases; they are
not hidden or converted to skips. The digest-pinned configuration executes all four
successfully with real copied raw inputs and unchanged assertions, not a mock
of source checks. The wrong digest fails at production binding validation;
the incomplete pair fails before any scientific test. Parent verification
recognizes the expected negative outcomes; its zero status must not be mistaken
for successful pytest results in those modes.

Each child blocks original-checkout, `/proc` and `/sys` canaries and forbids
subprocesses. The positive child loads 30 staged project origins with zero later
forbidden opens. All temporary copied staging is removed; private source
directories remain local. The 27 copied cases overlap the 200 local cases;
do not claim all 200 were run in the copied tree. Python-event guarding is
not OS containment, cross-host runtime proof or native biological re-admission.
[Machine-readable source and receipt pins](swiss_raw_pytest_handoff_20261002.json)
record the positive and negative evidence.

## Prior Platform Confirmation

Actual source-d03530ce run 36995416202 Linux job 110800898618 succeeds: **111
passes in 14.48s**, no skips/errors/failures. Its executed JUnit gate confirms
35 host-reader cases including the real live reader, plus the eight constructor
and 12 repeat-worker cases. The checkout SHA is verified in its actual log,
downloaded once at 10:27:55 UTC; artifact 11221322290 retains JUnit. This confirms
the preceding platform amendment, not this new pytest option patch or controlled
Threadripper timing. No local live host probe is run.

At 10:36:25 UTC the same run has successful Linux/docs/wheel; macOS 3.13/full
fail while 3.10/3.11/3.12 remain live. Inspect actual 3.13 job 110800898662 once
at 10:36:58 UTC: 14,024 passes/four failures/119 skips/zero errors/30 warnings
in 401.52s. All 34 portable reader cases pass; only its real Linux case skips,
with actual Linux execution verified separately. The four default raw exports
still fail because original paths are unavailable. Do not infer sibling counts,
new option confirmation or complete CI success.

Automatic public CI still has no raw data or these option values. The handoff
proves executable regression with legitimately supplied inputs, not automatic
public provisioning, data-rights clearance, new annotations/inference/scoring
or uncertainty estimation. Timing stays deferred without host contention polls,
questions, DGX access or unrelated job/service actions. Other-QfO uncertainty,
controlled resource evidence, rights, complete release and archival deposition
remain open; the original publication goal stays active and incomplete.
