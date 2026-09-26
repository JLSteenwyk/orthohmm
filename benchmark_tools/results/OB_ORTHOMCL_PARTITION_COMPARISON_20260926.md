# Retained OrthoMCL Run Differences

The [partition and family-count report](ob_orthomcl_partition_comparison_20260926.json)
compares the two retained native outputs without rerunning inference. File
identities match the previously recorded April and July hashes; duplicates
within groups, repeated gene assignments and repeated labels are rejected.
Labels and group ordering are ignored when comparing membership partitions.

April has 23,803 groups and July 23,804. Exactly 23,802 groups are identical.
The one April-only group contains ten genes. In July, seven of those genes form
one group and two form another; ENSP00000476467 is absent. Thus removing the
single April-only gene does not make the partitions equal. The report retains
the exact memberships of all three differing groups.

None of the ten genes in those changed groups appears in a retained RefOG.
Rescoring both partitions with the same 70 reference families and low-certainty
exclusions yields identical complete score objects, including all 70 family
TP/FP/FN/split/exact records:

| Statistic | April and July |
| --- | ---: |
| Weighted F1 (%) | 55.06534229215073 |
| Weighted precision (%) | 59.07805533272544 |
| Weighted recall (%) | 51.56306418084255 |
| Exact RefOGs | 12 |

The observed output differences lie outside the reference score's scope.
This establishes equality of this retained benchmark statistic, not general
partition equivalence, biological equivalence, shared historical commands,
input parity, or deterministic inference. The earlier summary's 216,950-gene
count matches April; the recent July readback contains 216,949. Both identities
remain explicit rather than silently exchanging their native coverage counts.

```bash
/usr/bin/python3 -S -m benchmark_tools.compare_ob_orthomcl_runs \
  --root . --output /tmp/ob-orthomcl-partition-comparison-new.json
python -m pytest -q tests/unit/test_compare_ob_orthomcl_runs.py \
  tests/unit/test_score_orthobench_partition.py
```

Ten comparator/scorer tests pass. No scores, inference defaults, runtime claims
or frozen publication configurations changed. Remaining provenance and
publication requirements are not completed by this comparison.
