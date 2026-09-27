# QfO Fresh Readback Failure 22330

## Terminal State

Slurm reports native job 22329 COMPLETED, exit 0:0, elapsed 01:53:08,
32 allocated CPUs. Readback 22330 FAILED, exit 1:0, after 16 seconds with
2 allocated CPUs. Neither job has been restarted. Native outputs remain
unadmitted for scientific comparison and cache reuse.

The retained log at
`benchmarks/work/qfo_fresh_phylogeny_20260927/readback-22330.log`
has SHA256 `e2548b4e99245059ac94529f1d4f54781cc0e95f1203758d685043f719b59483`.
It reports `ValueError: Changed complete input inventory` at line 77 of
`readback_qfo_fresh_phylogeny.py`, before scientific readback.

## Inventory Diagnosis

A direct record comparison of every planned FASTA against the input directory
found 78 expected files, no missing files and no changed records. The directory
contains one additional regular file, `staging_manifest.json`, 27,130 bytes,
SHA256 `07a890eb816944f946d46039559a6046c2a3664f033eaa0f30f35bb25b9c9ab8`.
The failing reader compares all regular files against a FASTA-only inventory.
The native worker passes the explicit planned FASTA filenames to the pipeline.
The staging manifest is not itself in the native plan's checked records; its
historical provenance must be investigated rather than asserted from this check.

## Native Receipt, Not Yet Scientific Validation

The native completion receipt has SHA256
`45e8f3164a0d3703302b632fdf368ffba0f6ce245401e8fff55f8358abadce53`.
It reports 351,739 candidate families, 24,264 reconciled, 327,475 bypassed,
366,068 root groups, 5,959,560 pairs, 26 species-tree marker families and zero
checkpoint reuse. Matching historical aggregate counts would not establish
equal predictions; the independent structural and pairwise checks remain undone.

GNU time reports 1:53:07 wall, 192,053.37 user seconds, 9,725.70 system seconds,
and 7,661,580 KiB maximum process RSS. These are shared-host descriptive
measurements, not dedicated-machine timing or peak process-tree memory.

## Next Actions

Preserve the failed audit plan, submission and log. Correct and test the reader's
inventory policy with explicit metadata handling, without silently ignoring
arbitrary added files. The canonical reader contains the same inventory pattern
and needs the same correction. Freeze a new audit plan and submit a distinct
readback job only after testing; update the downstream gate to require that
successful audit while retaining the failed attempt. Do not edit the original
native plan, rerun native inference, transfer historical scores, or submit the
canonical run before admission succeeds. Publication readiness remains unproven.
