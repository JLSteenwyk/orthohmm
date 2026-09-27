# FastOMA OrthoBench Provenance Readback

The [retained-file audit](ob_fastoma_provenance_20260926.json) connects the
FastOMA OrthoBench comparison to its staged proteins, supplied tree, final
collection task and published output, retaining successful and failed task
records. No native inference or scoring was rerun.

## Inputs, Tree and Output

All 12 staged proteomes match the frozen OrthoBench gene IDs and sequences
exactly, totaling 251,378 proteins. Matching uses unique species-level gene-ID
sets followed by exact sequence comparison, not filename or file-order
assumptions. The audit records the original-to-short-species-name mapping.

The input species tree and checked output have the same rooted clades and
all 12 species. The checked tree adds branch lengths, so this is topology
agreement, not byte identity. The supplied tree's original provenance is not
established here, and this was not independent FastOMA species-tree inference.

| Output | Groups | Assigned genes | Role |
|---|---:|---:|---|
| OrthologousGroups.tsv | 18,972 | 130,010 | Scored final groups |
| RootHOGs.tsv | 19,172 | 197,878 | Separate diagnostic, not substituted |

The published final table is byte-identical to the successful collection
task's table already bound to the retained score. Its SHA-256 is
`7599b5aac0c871be6114e2f845b19c18f736839a655a1aa535a19539918f1dbd`.
Group-table membership is unique, nonempty and restricted to the input genes.
Higher root-HOG coverage must not be used as final-group coverage.

## Retained Execution History

The log reports FastOMA pipeline version 0.3.5 and a Docker-profile Nextflow
command. The exact command and all task command/wrapper/log hashes are retained.
There are **100 task exit records: 90 successful and 10 with exit code 137**,
matching the workflow summary. The failed attempts are preserved. Exit 137
alone does not establish out-of-memory as the cause, and the completed workflow
must not be described as a failure-free run.

The workflow reports 50m 23s duration and 27.3 CPU hours (0.2% failed). These
rounded workflow-reported observations are not externally verified controlled
timing. Peak RSS is unavailable in this audit. No speedup or memory-efficiency
claim follows. Historical container/database identity and full input-consumption
provenance remain unproven despite current file consistency.

## Verification

Twenty-two focused FastOMA/OrthoFinder provenance tests pass, covering input
identity, mismatched residues, ambiguous species mapping, native output types,
duplicate/unknown genes, topology changes and preservation of failed tasks.
All evidence hashes are checked again after readback.

```bash
python -m benchmark_tools.audit_ob_fastoma_provenance \
  --root . --output /tmp/ob-fastoma-provenance-new.json
```

This supplements the retained FastOMA score without changing it. It does not
complete the remaining cross-method provenance, controlled resource or release
requirements, and does not establish publication readiness.
