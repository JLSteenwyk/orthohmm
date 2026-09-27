# QfO OrthoFinder Provenance Consolidation

The [machine-readable receipt](qfo_orthofinder_provenance_consolidated_20260926.json)
joins both current QfO OrthoFinder rows to their shared, hash-pinned native
admission, execution, command, timing, inputs and pair conversions. It rechecks
185 records, including both submitted/filtered pair files and 78 copied FASTAs.
It does not rerun all 465,149 native-admission record checks or QfO scoring.
Internal sequence parity is inherited from the pinned prior admission, which
reports 984,137 identical sequences across 78 species and no differences.

Both rows use OrthoFinder 3.1.5, native job 21706, on bizon with 32 allocated
CPUs. The command was `orthofinder -f INPUT -t 32 -a 32 -S diamond`, with its
absolute executable and input paths retained in the receipt. Both the native
log and GNU-time command agree with the recorded execution command.

| Output | Submitted and retained pairs | Pair-conversion wall seconds | Separate inference time |
| --- | ---: | ---: | --- |
| Full native phylogenetic ortholog pairs | 14,215,382 | 193.521242 | 42,018 s full-run GNU time |
| Pre-phylogeny MCL group-derived clique pairs | 163,277,439 | 249.173882 | Unknown; extracted from the same full run |

The checkpoint has 167,209 groups. Reference filtering removed zero pairs and
preserved file hashes for each output. Pair volume is not protein coverage and
does not itself measure accuracy. No score changed: the six QfO endpoints and
secondary means are copied exactly from the frozen eight-method manifest.

Full-run GNU time reports 453,857.43 user seconds, 9,128.89 system seconds and
6,479,368 KiB maximum process RSS. Scheduler elapsed time is 11:41:56, whereas
the timed native command is 11:40:18; these scopes must not be interchanged.
Conversion intervals come from the converter's epoch timestamps. Maximum
process RSS is not simultaneous aggregate process-tree memory. Shared-host
measurements support no controlled comparative speed or memory ranking.

The previously admitted native runtime is linked through the frozen command
plan. This consolidation does not repeat its transitive executable audit.
Other tools' QfO provenance, controlled scaling, and publication requirements
remain separate work. Neither historical results nor production defaults changed.

Reproduce to a new destination:

```bash
python -m benchmark_tools.consolidate_qfo_orthofinder_provenance \
  --repo . --output /tmp/qfo-orthofinder-provenance.json
```

Thirty focused provenance/binding and shared resource-parser tests pass.
