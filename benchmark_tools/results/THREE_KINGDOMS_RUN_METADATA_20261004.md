# Retained Three Kingdoms Command And Resource Supplement

Source `3ef56766` executes once successfully using the private stdlib controller.
The [actual supplement](three_kingdoms_run_metadata_20261004.json) recovers seven
historical rows, six GNU-time command/resource records and two OrthoHMM metrics
records. The contemporary matched Sonic row is already bound in the
[selected all-tool register](all_benchmark_provenance_20261004_v2/register.json);
its older historical directory is deliberately excluded. No inference,
conversion, scoring or timing replicate is rerun.

```bash
python -B -m benchmark_tools.audit_three_kingdoms_run_metadata \
  --repo . \
  --output benchmark_tools/results/three_kingdoms_run_metadata_20261004.json
```

The helper verifies each selected retained prediction's actual bytes/hash
before joining its result-directory metadata. It checks50 direct records
plus source and five helper identities. Binding is a retained directory
association with a prediction identity, not immutable historical execution
attestation. Version text/executable names are recorded, not universally
attested runtime identities. Raw input/source manifests are inherited.

## Recovered Measurements

| Method/stage | Seconds | CPU user/system seconds | Memory scope |
| --- | ---: | --- | --- |
| OrthoHMM high sensitivity metrics | 6556.426547 | 189932.597022 / 92.579848 | Sampled tree RSS,11941715968 bytes |
| OrthoHMM phylogeny metrics | 8222.745564 | 231374.661843 / 7087.538037 | Sampled tree RSS,11784400896 bytes |
| OrthoFinder full GNU time | 6332 | 69581.08 / 1890.13 | Maximum process RSS,4178948 KiB |
| OrthoFinder MCL checkpoint conversion | 0.84 | 0.68 / 0.16 | Conversion maximum process RSS,150084 KiB |
| ProteinOrtho GNU time | 568.44 | 8847.27 / 55.67 | Maximum process RSS,2698752 KiB |
| FastOMA GNU time | 6522 | 84.99 / 99.46 | Driver maximum RSS,621124 KiB; not aggregate Docker-task memory |
| OrthoMCL BLAST-to-BPO conversion | 547.24 | 534.73 / 12.45 | Conversion maximum process RSS,72392 KiB |
| OrthoMCL recovered mode4 | 2276.32 | 6127.68 / 3109.58 | Downstream maximum process RSS,16676844 KiB |

Keep the previous8225/569/6523-second wrapper observations separate from
the8222.745564/568.44/6522-second native/metrics intervals. They are not new
replicates or replacements. OrthoFinder's derived4884-second checkpoint
timestamp is not its0.84-second conversion nor separately measured inference.
OrthoMCL's accepted integer stage sum remains
138194 BLAST +548 conversion +175 indexing +2276 downstream =141193 seconds.
Do not substitute downstream2276.32 seconds for that stage sum or full search.

GNU-time CPU covers the child and waited-for descendants, not necessarily
concurrent container tasks. OrthoHMM CPU is retained collector output without
fresh process-tree/cgroup validation. Distinct memory scopes cannot be ranked
as comparable peaks. These are shared-host historical observations, not
isolated performance or causal efficiency estimates.

## Chronology And Scope

Ordered metadata retains recovery jobs20900,20901,20902 and20903. The earlier
attempts are not silently overwritten by a last-value dictionary. Their exit0
records do not override the prior scientific invalidations. The accepted
downstream scope remains the final recovered mode4 run. The empty OrthoMCL
`time.log` is explicitly recorded as empty, not a zero-duration full pipeline.

Launch-source files record599dda43 for high sensitivity and6bb08e5b for
phylogeny; harness-at-completion records63aeb8a8 and9c293177 respectively.
Both harness snapshots mark the worktree dirty. Review of historical
`benchmark_production.py` at63aeb8a8 confirms Git/source/input manifests are
collected after the child finishes. These fields have different temporal
scope; do not claim the completion commit uniquely identifies launch code or
that post-run input hashes prove historical consumption.

The selected phylogeny group identity matches the harness's retained clustered
partition entry, not the separately recorded final singleton-added group
entry. Preserve the existing BUSCO group co-membership endpoint and do not
reinterpret it as native pairwise ortholog scoring. No endpoint or score changes.

## Validation And Identities

Twenty new focused cases plus the existing30 register cases pass:50 total,
1.32s pytest display,1.279s JUnit suite time, no errors/failures/skips.
Tests cover repeated recovery keys, malformed metadata, failed/unbound metrics,
invalid resource numbers/argv, distinct harness/source scopes, replaced Sonic
exclusion, changed evidence, failed native time and existing-output preservation.

Independent stdlib readback rechecks all50 direct records and six source/helper
pins; compares metrics commands, measurements and inherited manifests with the
original JSON; verifies the selected output bindings, ordered recovery IDs,
exact stage arithmetic and fractional GNU-time intervals. Readback passes.

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| Actual supplement | 83910 | `bb6d910021583da61ba4892690c5bf50796fcf81cafcfaeb0d09742e41b8a093` |
| `audit_three_kingdoms_run_metadata.py` | 9456 | `669c4ab84acfd6fe1df73e37940780e6abc09be3fadc64bb399af855e6b9d64c` |
| `tests/unit/test_audit_three_kingdoms_run_metadata.py` | 6575 | `0348d4d4c019fca45d10d30f6cedaca23b68355997e3c24c729b97c0a32ff9fd` |
| Local `benchmarks/work/three_kingdoms_run_metadata_tests_20261004.xml` | 7090 | `1863141c348dde004c827e41b0f6027e37886fcffdcbc076f2c0e4ae6ef7ef07` |

No universal historical runtime/input attestation or complete publication goal
is certified. Unresolved scientific and transitive delivery requirements in
the full-goal audit remain. Preserve existing register/manuscript/audit/rc1
bytes and use this newer supplement alongside them. No owned native job is live.
