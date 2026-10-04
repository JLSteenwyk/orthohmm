# Actual All-Tool, Three-Dataset Register

Execute the corrected join from source `b13883bacc1291a1a825e3567987857bc62fbd38`.
The selected [register](all_benchmark_provenance_20261004_v2/register.json)
contains exactly24 method/dataset rows: eight methods each on OrthoBench,
corrected QfO and Three Kingdoms. Its [resource table](all_benchmark_provenance_20261004_v2/resources.tsv)
has29 entries because available inference and conversion costs are separate.
The [compact index](all_benchmark_provenance_20261004_v2/register.md) links the
resource scopes and retained input notes. This is consolidation, not new
inference, scoring, raw-input admission or historical execution attestation.

## Execution And Readback

The private stdlib controller runs with Python/site/preload environment
overrides removed, bytecode disabled and a new absent cache-prefix path:

```bash
python -B -m benchmark_tools.assemble_benchmark_provenance \
  --repo . \
  --output benchmark_tools/results/all_benchmark_provenance_20261004_v2
```

Actual exit0 reports `all_benchmark_provenance_register_partial`,24 rows
and32 directly checked metadata/command/time records. The source and four
imported helper identities are also recorded. Output/input identities copied
from existing admissions are explicitly inherited, not newly raw-file hashed.

Independent stdlib readback, without importing the assembler, verifies the
complete Cartesian method/dataset inventory, every individual score and QfO
secondary mean against the selected current score manifest, and all prediction
bindings against their score/conversion records. It rechecks32 direct-record
and five source/helper byte counts and SHA256 identities. All29 TSV entries
match the structured register including scope, units and null values.
Additional checks confirm exactly one conversion entry per OrthoFinder QfO
row,78 actual Sonic QfO input records rather than relation-table names, the
current Three Kingdoms Sonic4533-second run, and the non-full OrthoBench cached
319.467197-second and checkpoint0.42-second stages. Readback passes.

The retained final JUnit receipt reports30 tests, no errors/failures/skips,
1.019s suite time (pytest displayed1.04s). These are focused synthetic adapter,
binding, refusal, fresh-output and resource-scope tests, not additional native
benchmark replicates. No owned Slurm job is live at this checkpoint.

## Resource Interpretation

| Dataset and method | Seconds | Bound scope |
| --- | ---: | --- |
| QfO OrthoFinder full | 42018 | Full native inference; same run also supplies the MCL checkpoint |
| QfO SonicParanoid | 11423 | Full native inference |
| QfO ProteinOrtho | 4709 | Full native inference |
| QfO FastOMA | 14094 | Native workflow wall time; GNU-time CPU/RSS concern the driver and waited descendants, not aggregate Docker tasks |
| QfO OrthoMCL | 7455 | Recovered mode4 downstream only; excludes BLAST/BPO |
| Three Kingdoms SonicParanoid | 4533 | Contemporary matched-input native run, not historical4479-second run |

Do not rank efficiency from these values. Allocations, input histories and
memory scopes differ. GNU-time maximum-process RSS, sampled process-tree sums,
driver RSS and stage-sum observations retain their distinct labels and units.
Three Kingdoms OrthoMCL's sampled sum double-counts shared fork pages. QfO
OrthoHMM full search-to-current-output costs remain unknown; incremental
reconciliation timing is not a substitute. The replacement shared-host
Threadripper panel remains separate and is not assigned to these score runs.
Contention effects remain unknown and potentially method-dependent.

## Development History

Preserve the first actual invocation's pre-output failure: five OrthoBench
prediction pins used the supplemental score schema. Add missing/conflicting
pin checks and schema coverage before the27-case successful first export.
Readback of that unselected `all_benchmark_provenance_20261004` export found
duplicate OrthoFinder conversion entries and a FASTA-suffix heuristic that
misclassified Sonic relation tables as inputs. The corrected source uses
explicit native input manifests, emits conversion once, and retains known
cached/checkpoint costs and proper memory measurement scopes. Three additional
cases bring the final test count to30. Retain older outputs and local JUnit
receipts as history; do not overwrite or select the defective first export.
These were reporting/schema defects, not native/scoring changes.

## Selected Identities

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| `assemble_benchmark_provenance.py` | 19979 | `7e1b8c36ae61277d54f42808855e80738dfe2225c3fc52381b354bece11675bc` |
| `tests/unit/test_assemble_benchmark_provenance.py` | 12258 | `f1df86dde1c6411eb9a2aa58e34f7832f51439b33dd3b6110ba3927b1a6e0bbb` |
| Local `benchmarks/work/all_benchmark_provenance_tests_v4_20261004.xml` | 4304 | `06d3e386c2008ffac6a966e19615cb880cc1e504232d0ef5050fea041fdca439` |
| Selected `register.json` | 458287 | `48be4b216f9f3776354246894c8c2a38f790fd4f075528e88f56b230a0a3c47c` |
| Selected `register.md` | 4502 | `de63981db48168aad114ddcaf3b2d440604c1e7f5c3c669608b19b2c313dfeb8` |
| Selected `resources.tsv` | 3709 | `674a9bb16e749cb57ed3bde44931a443d91dd2f3b861d4548ffa93698e9486b4` |

The nine source receipt identities and all32 directly checked records are in
the structured register. Preserve prior immutable audit/manuscript/rc1 bytes;
this newer metadata result does not retroactively change their evidence pins.

## Remaining Work

This advances goal1.2/1.4 but does not certify complete transitive provenance,
universal runtime version attestation or historical input consumption. Recover
available historical commands/runtime records only from bound retained
evidence, leaving missing fields explicit. Unresolved QfO uncertainty, partial
error-stratum/mechanism validation, development-family inventory, per-ablation
full-pipeline costs and transitive study delivery remain scientific/delivery
gaps. The full publication goal stays active and incomplete. No quiet-window,
DGX, new family-disjoint or host-certification gate is introduced.
