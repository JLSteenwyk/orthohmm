# SonicParanoid Native Table Ownership Audit

Validated the original-release native tables under
`qfo_benchmark/results/sonicparanoid/run/sp_out/runs/sp2_25426142059_default_32cpus_ml075/species_to_species_orthologs`
against all canonical FASTAs in `qfo_benchmark/input`.

| Quantity | Verified count |
| --- | ---: |
| Input proteomes | 78 |
| Unique input accessions | 976,504 |
| Species-pair tables | 3,003 |
| Native table rows | 4,994,833 |
| Raw relations | 15,028,392 |
| Distinct converted pairs | 15,022,677 |
| Repeated relations removed by existing converter | 5,715 |

Every expected input-species pair has exactly one table. Header, row counts,
member-list structure, accession ownership, accession-alias uniqueness
within a member list, and finite confidence values pass. Empty tables are
permitted. Repeated relations across rows retain the existing converter's
deduplication behavior; raw = distinct + repeated is checked explicitly.

The audit recorded SHA-256 and size for all 3,003 tables and 78 input files
before traversal and rechecked all 3,081 records afterward. The table-file
inventory was also unchanged. This run returned successfully; these checks
do not establish historical workflow completeness or biological accuracy.
No converter, original prediction or score was changed. Corrected outputs
still require their own successful execution binding and graph/table audit.

## Reproduction and Tests

```python
from pathlib import Path
from benchmark_tools.validate_sonicparanoid_tables import validate_tables

directory = Path("qfo_benchmark/results/sonicparanoid/run/sp_out/runs/sp2_25426142059_default_32cpus_ml075/species_to_species_orthologs")
print(validate_tables(directory, sorted(Path("qfo_benchmark/input").glob("*.fasta"))))
```

Sixteen focused tests passed across `test_validate_sonicparanoid_tables.py`
and the existing `test_sonicparanoid_to_pairwise.py`. They cover member
ownership, unknown accessions, nonfinite confidence, count mismatches,
duplicate aliases, repeated relations, empty tables, an entirely omitted
species, and reversed duplicate species tables.

## Progress Ledger

Previous turn: progress, manuscript/input-release claim synchronization
pushed as `3b59699`. Current turn: implemented and tested full-input table
validation and successfully audited retained native tables.

QfO assessment `21711` remained running at 42:01 with admission `21712`
pending; no score was inferred from runtime or logs. Corrected inference
and the dedicated timing array remain in progress. Next: bind corrected
SonicParanoid tables to their terminal execution record, validate them,
then perform separately frozen conversion/filtering/scoring. The complete
publication goal remains open.
