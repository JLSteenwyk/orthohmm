# Complete BPO Content Consistency Check

`audit_orthomcl_bpo_content.py` streams the source m8 table and converted BPO
together. It checks each retained directed pair's similarity ID, sequence
identities and lengths, first-HSP E-value, subject-span-weighted truncated
percent identity, and every HSP span. It detects omitted/extra/reordered BPO
records, noncontiguous source blocks, malformed source fields, inconsistent
alignments and changed artifacts. Input and helper hashes are rechecked before
writing a new report; existing reports are never overwritten.

Cutoff-excluded pairs are counted while their source HSPs are still read and
validated. The nominal cutoff remains `1e-5`. First-HSP significance is
intentional: the retained BioPerl `Bio/SearchIO/blasttable.pm` sets
`Hit_signif` when starting a query/hit block (reviewed at lines 254-270), and
native OrthoMCL uses that hit significance for filtering. The checker does not
silently replace this rule with the minimum E-value across HSPs.

The sequence-length lookup is shared with the converter; the expected BPO
record construction and source-to-output comparison are separate. Identity
arithmetic follows the reviewed native-compatible binary floating convention,
including integer truncation. Therefore this is full-file consistency
verification, not proof against every shared algorithmic error. The native
OrthoMCL/BioPerl parity probe remains complementary evidence.

## Retained Evidence

Both the native and streaming BPOs from the retained parity fixture passed:

| Quantity | Count |
| --- | ---: |
| Input proteins | 4 |
| Source HSP rows | 8 |
| Source directed-pair blocks | 7 |
| Cutoff-excluded blocks | 1 |
| Verified BPO pair records | 6 |

The saved reports are `orthomcl_native_bpo_content_probe_20260918.json` and
`orthomcl_streaming_bpo_content_probe_20260918.json`. They bind the existing
native/streaming outputs under
`benchmarks/work/orthomcl_bpo_parity_probe_v2_20260918`; no fixture or historical
production output was rewritten.

All 48 focused content/index/converter tests passed with installed native
parity testing enabled. Tests alter every BPO field, second-HSP spans,
cardinality/order, source block contiguity and input hashes. They also verify
the first-HSP cutoff rule when a later HSP has a lower E-value.

## Remaining Integration

No corrected production BPO exists yet. The converter, complete-content check
and native index validator must be placed behind the queued search admission,
with frozen Perl/BioPerl/runtime provenance and terminal execution checks.
Only then can native clustering and final-group validation proceed.

Memory retains sequence lengths, query IDs and one pair's HSP spans rather
than the complete hit table/BPO. Index validation has separate memory costs.
This audit does not itself establish search provenance, query success,
clustering correctness, benchmark accuracy or publication readiness.
