# Historical FastOMA OrthoXML Audit

The complete historical QfO FastOMA OrthoXML, root-HOG table and native pair
file passed identifier and membership consistency checks against the frozen
original-release input FASTAs. This is not the corrected-release run.

| Quantity | Count |
| --- | ---: |
| Input species | 78 |
| Input proteins | 976,504 |
| Proteins declared in OrthoXML | 967,184 |
| Input proteins not declared in OrthoXML | 9,320 |
| Proteins referenced in HOGs | 588,873 |
| Declared proteins not referenced in HOGs | 378,311 |
| Root HOGs | 55,486 |
| Proteins represented in native pairs | 584,445 |
| HOG members without a native pair | 4,428 |
| Native pair rows checked | 15,320,615 |

All species declarations match the input owners and XML taxonomy. Gene IDs and
accessions are unambiguous; gene references are defined and occur once. Group
identifiers and taxon references are valid, numerical scores are finite, and
no root HOG is empty. `RootHOGs.tsv` exactly matches all XML root memberships.
Every native pair consists of known cross-species proteins in the same root
HOG. The separate distinct-pair audit established that all these pairs are
unique; this audit does not repeat that global uniqueness assertion.

The parser uses streaming ElementTree events and releases parsed XML elements.
It retains identifier and membership maps, not the complete XML document.
Root-HOG names follow the native `collect_subhogs.py` prefix convention before
the underscore. No orthology relations are inferred from root co-membership:
native pairs remain the evaluated predictions.

Evidence: `fastoma_orthoxml_audit_20260918.json`, SHA-256
`6a6124010346eaffd0bf0169dbbabe9a713c8398066e4757fd3859649b822fe7`.
The report binds the complete inputs and checked outputs by checksum; files
were rechecked after analysis. No historical predictions or scores changed.

Reproduce with a fresh report path:

```bash
python benchmark_tools/audit_fastoma_orthoxml.py --root . \
  --output benchmarks/work/fastoma_orthoxml_audit_fresh.json
```

Validation: 68 focused XML/native-pair tests passed, followed by the full real
historical-data audit. Tests include malformed hierarchy, foreign/duplicate
identifiers, unresolved/repeated references, mismatched taxonomy, nonfinite
scores, root-table differences and native pairs outside root scope.

## Interpretation Limits

These are native-output coverage counts, not reference-relative recall. The
checks do not identify why proteins were omitted or prove that omissions are
biologically appropriate. Root co-membership is necessary but not sufficient
for orthology. Complete Nextflow task-chain admission and biological accuracy
remain separate requirements. The reusable validators will support corrected
FastOMA native-output admission after that inference job finishes.
