# OrthoFinder OrthoBench Provenance Readback

The [machine-readable audit](ob_orthofinder_provenance_20260926.json) connects
the retained full-pipeline and sequence-only OrthoBench prediction files to
the available input, native checkpoint, version/command and resource records.
All seven checks pass. This is retrospective retained-file consistency, not
independently attested historical execution or a reconstructed binary history.

## Inputs and Predictions

- All 12 staged FASTA files match the hash-pinned frozen OrthoBench input
  inventory byte-for-byte.
- All 251,378 internal processed protein sequences match those FASTAs exactly
  after applying `SequenceIDs.txt`; species mappings and IDs also match.
- The retained full output has 35,454 groups covering every input gene once.
  Its hash matches the full prediction used by the publication comparison.
- The retained sequence-only output has 34,230 groups covering every input
  gene once. Translating the native MCL checkpoint with its ID map reproduces
  the complete partition, independent of row order and group labels.
- The sequence-only prediction hash matches the previously admitted comparator
  readback. No benchmark scores or predictions were changed.

The checkpoint check reuses the existing converter's MCL parser, rather than
claiming independent parser validation. Full-output semantics are supported
by the separately recorded [version-specific source audit](ORTHOBENCH_UNCERTAINTY_PROTOCOL_20260916.md),
not inferred from the name `Orthogroups.txt`.

## Command and Resources

The native log reports OrthoFinder **3.1.5**, successful completion, and the
same command as the GNU-time log:

```text
orthofinder -f <retained OrthoBench input> -t 32 -a 8 -S diamond
```

Exact executable and input paths are retained in the JSON. The conversion
command names the same native checkpoint, mapping and sequence-only output.

| Recorded stage | Elapsed seconds | User CPU seconds | System CPU seconds | Maximum process RSS (KiB) |
|---|---:|---:|---:|---:|
| Full native inference | 3421.07 | 64429.47 | 3743.77 | 4787980 |
| Sequence-checkpoint conversion | 0.42 | 0.34 | 0.07 | 85744 |

Both time logs report exit status zero. **Sequence-only inference runtime is
unknown**, not 0.42 seconds: the comparator was extracted from a checkpoint
of the full run, not independently timed through inference. These historical
process measurements do not establish controlled efficiency, simultaneous
process-tree memory, or complete preprocessing/conversion/scoring totals.

## Verification

The audit records hashes before and after reading every input/evidence file.
Twenty-six focused provenance, converter and resource tests pass, covering
changed sequences, missing/extra/duplicated IDs, mapping errors, exact partition
coverage, command parsing and duplicate group labels. The retained evidence
passed without inference or scoring reruns.

```bash
python -m benchmark_tools.audit_ob_orthofinder_provenance \
  --root . --output /tmp/ob-orthofinder-provenance-new.json
```

This supplements two OrthoBench rows; it does not close the remaining tools'
native provenance or the dedicated-machine resource requirement. Current
binary/dependency hashes cannot retroactively prove which binaries executed
in July. Publication readiness remains unproven.
