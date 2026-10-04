# Complete Requirement Inventory, Incomplete Publication Goal

The [human-reviewed register](PUBLICATION_REQUIREMENT_AUDIT_20261004.md)
and [machine-readable audit](publication_requirement_audit_20261004/audit.json)
cover all47 authoritative requirements: five scope paragraphs,36 numbered
bullets, five engineering bullets and one wrapped completion paragraph.
Thirty-seven have bounded support, nine are partial, and the full completion
criterion is unmet. This is not a percentage-complete estimate or a new
scientific admission. No goal-complete or publication-ready status is set.

The stdlib generator enforces exact requirement coverage and pins55 direct
evidence files, not their complete transitive scientific dependencies. An
independent line-based readback checks the36 numbered/five engineering bullets,
five scope IDs, exact joined completion paragraph and all55 evidence identities.
No same-parser-only completeness claim. Eighteen focused cases pass in0.27s
(JUnit0.246s), including adjacent bullets/headings, wrapped text, missing/extra
requirements, duplicates, invalid status/evidence and preservation of outputs.
Even a wholly `supported` synthetic register cannot emit goal certification.

The first actual generator invocation refuses the real goal's adjacent
completion heading; the second refuses its incorrectly merged bullet inventory.
Neither writes an audit. Both formatting cases are added to tests before the
successful real47-row invocation. Do not alter the authoritative goal's layout
to make the parser pass. Preserve the earlier16/17-case local receipts as
development history; the18-case receipt is the final tested source.

```bash
python -B -m benchmark_tools.audit_publication_requirements \
  --goal /absolute/path/to/authoritative-goal.txt \
  --review benchmark_tools/results/PUBLICATION_REQUIREMENT_AUDIT_20261004.md \
  --output /absolute/new/requirement-audit
```

Use the original goal or the exact retained `goal.txt` snapshot. This performs
inventory/identity checks on a human assessment; it does not infer scientific
sufficiency from hashes, green tests or existence of receipts.

| Current artifact | Bytes | SHA256 |
| --- | ---: | --- |
| `audit_publication_requirements.py` | 5781 | `9b8a0e46d884a1d9b831666328d793ba839b58c7733d8aa460ac0364d14be60f` |
| `tests/unit/test_audit_publication_requirements.py` | 5742 | `127675425efe76994e0807f45b494b87ded77ba61f0cf8f54cb08d26ae5c122e` |
| Human register | 22482 | `fca9cf40c4223ac3e43db649f91030664fc26b3f60b3652fe0f8005bf6ce393c` |
| `publication_requirement_audit_20261004/audit.json` | 56294 | `dd2ce3aa7849995f3c9780c570063a3cd2476964410a0e1ffa51748e8705f9e8` |
| `publication_requirement_audit_20261004/goal.txt` | 14267 | `4028ccc3863abae74c26973cf577fa0ffcc717909d326f2811a5c4bf02e237b7` |
| Local `benchmarks/work/publication_requirement_audit_tests_adjacent_20261004.xml` | 2652 | `a09d22507f6385ee7e29e37f6f82732a25cfe068a5c9e6328082c3429e9a2c90` |

## Next Scientific And Delivery Work

Consolidate the all-tool/all-dataset provenance and resource gaps from existing
receipts first (1.2/1.4), without claiming historical input consumption from a
present-day checksum. Nine partial requirements are1.2,1.4,2.1,3.4,4.1,4.3,4.4,
7.4 and7.5. Retain all missing/incremental per-ablation cost labels, unresolved
native QfO intervals, incomplete development-family inventory and descriptive
feature/mechanism limitations. Only independently justified additional analyses
can close those gaps; an arbitrary resampling unit or replay cost cannot.

The local rc1 candidate and32-page presentation are real completed deliveries;
all-tool transitive reproduction and final release metadata/format are not.
Supplement the candidate with the current audit rather than rewriting its
historical payload bytes. Public deposition/submission remain explicitly
unexecuted, not fabricated identifiers or newly imposed permission gates.
No quiet window, DGX, OS/security certification or new family-disjoint gate is
required. No owned native job is live and no completed successful scientific
analysis or27-identity panel attempt was restarted for this audit.
