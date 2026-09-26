# OrthoBench Comparator Prediction Readback

The [readback report](retained_ob_comparator_readback_20260926.json) supplies
current prediction-file identities for the five comparator rows whose retained
summary objects omit direct hashes. All five reproduce retained precision,
recall and F1 within 1e-8 percentage points and exact-RefOG counts exactly.
Seventy reference families and their retained low-certainty exclusions were
read under the existing official-formula implementation. Ninety-two evidence
records were checked before/after their use, including helper files.

| Method | Readback F1 (%) | Input |
| --- | ---: | --- |
| OrthoFinder sequence checkpoint | 58.705963184 | `orthofinder_v3_sequence_only_20260830/orthogroups.txt` |
| SonicParanoid | 46.757609120 | `orthogroups_sonicparanoid.txt` |
| ProteinOrtho | 45.057346674 | `orthogroups_proteinortho.txt` |
| FastOMA final OGs | 30.906942172 | native two-column `OrthologousGroups.tsv` |
| OrthoMCL | 55.065342292 | retained July 25 `all_orthomcl.out` |

The initial generic `orthogroups_fastoma.txt` candidate instead gives
50.818044059% F1 and 21 exact groups, matching the documented root-HOG
diagnostic rather than the final-OG endpoint. It was not substituted into the
comparison. The final two-column table is parsed by the existing FastOMA
adapter; OrthoMCL uses the existing annotated-member adapter. The report
retains exact absolute paths, hashes, parsers and per-family counts.

```bash
/usr/bin/python3 -S -m benchmark_tools.audit_retained_ob_comparators \
  --root . --output /tmp/retained-ob-comparator-readback-new.json
python -m pytest -q tests/unit/test_audit_retained_ob_comparators.py \
  tests/unit/test_score_orthobench_partition.py \
  tests/unit/test_normalize_three_kingdoms_orthogroups.py
```

Thirteen focused comparison/scorer/adapter tests pass. No native inference,
historical conversion or official evaluator executable was rerun here. Matching
the statistics establishes current-file readback equivalence, not proof that
these exact files were consumed historically. Commands, versions, native
conversion provenance, resource evidence and historical input consumption still
need separate consolidation. Existing score tables are unchanged; this report
supplements rather than retroactively rewrites their provenance.

The subsequent [upstream-function cross-check](RETAINED_OB_UPSTREAM_CROSSCHECK_20260926.md)
returns exactly identical precision, recall and F1 for all five rows. It also
documents a one-gene inventory difference between the retained April and July
OrthoMCL outputs; current readback uses July and is not a partition-identity claim.
