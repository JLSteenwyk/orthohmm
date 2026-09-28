# TreeFam-A 7-Derived Subtrees Recovered

The public [Selectome download page](https://selectome.org/download) links an
[archive index](https://selectome.org/ftp/MySQL/) containing
`selectome_03-TF7A__mysql5.0.sql.zip`. Downloaded this archive, the index,
[published checksums](https://selectome.org/ftp/MySQL/SHA256.sum) and the FTP
parent index into `benchmarks/work/treefam_selectome_search_20260928/`.
No maintainer or other person was contacted.

The 2,429,072-byte archive has SHA256
`4047be74d730e7057633b5282abfe999cc0202dd42eb194b8efc69586c2745d1`,
matching the publisher's checksum. Its sole member is the 15,917,948-byte SQL
dump, SHA256
`13785609901d4fe0d0f04eb3edecd169e2920772a7f09c6a465d2d7e209aba5e`.
All ZIP CRC checks pass. The embedded news table explicitly identifies a
TreeFam-A 7-based release, dated May 18, 2009; the dump footer is June 2, 2009.
These statements establish the archive's declared origin, not identity to QfO's
original reference-generation files.

## Inspected Contents

The [machine-readable inventory](selectome_tf7a_inventory_20260928.json) retains
all seven table names, the release news, every subtree's family/taxon/number,
NHX hash, species set, leaf count and exact duplicate-label counts.

- 1,211 parsed NHX subtrees representing 1,063 family accessions.
- All subtree taxa are Euteleostomi; leaves span 13 species labels.
- 20,781 leaf occurrences, not necessarily distinct genes across subtrees.
- No exact duplicate leaf labels within a subtree; eight subtrees contain
  labels differing only in case. The initial taxon-normalizing parser rejected
  the first such collision. The final inventory preserves leaf labels as node
  labels, without case folding or merging them.
- No `treefam2reference` filename string is present in the full SQL payload.
  Tables are Selectome branch, codeml, news, position, subtree and summary tables,
  plus taxonomy links; no original QfO mapping table was recovered.

SQLGlot 30.20.0 parsed the selected MySQL INSERT statements as literal data;
DendroPy parsed the NHX strings. The dump was never imported or executed as SQL.
SQLGlot was installed only into the local search directory, not a frozen
scientific runtime. Seven parser tests pass, including escaped strings and
rejection of functions, additional statements, wrong tables and nonliteral
inserts. A second full parse reproduces every inventory field.

## Interpretation And Reproduction

This corrects the earlier search boundary: the first Selectome paper describes
TreeFam 4/6, but its later public archive does preserve TreeFam-A 7-derived
subtrees. It does not recover an authenticated complete TreeFam-A release-7
collection or the QfO mapping. Vertebrate subtrees cannot silently replace full
families. No mapping suitability, pooled-reference reconstruction, benchmark
family assignments or family-level uncertainty has been established.

The executable inspector pins the archive checksum and expected member size.
With SQLGlot 30.20.0 and DendroPy available, run:

```bash
python -B -m benchmark_tools.inspect_selectome_tf7a --directory benchmarks/work/treefam_selectome_search_20260928 --output /tmp/selectome_tf7a_inventory.json
```

Raw downloads remain outside Git. The public index also lists a TF7AB-derived
dump; it is not inspected or admitted here. No original mapping, newer release,
inferred family labels or derived subtrees were substituted into the benchmark.
