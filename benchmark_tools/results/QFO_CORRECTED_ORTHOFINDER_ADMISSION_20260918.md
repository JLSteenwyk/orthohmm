# Corrected OrthoFinder Native Admission

`admit_qfo_corrected_orthofinder.py` is the independent native evidence gate
for corrected primary task `21706_1`. It does not convert predictions or score
accuracy, and its existence does not mean corrected OrthoFinder has finished.

The gate requires terminal scheduler success, 32 CPUs on bizon, the pinned
primary plan SHA-256 `dbebd3a6915fddeb2b89e0e591a2ac5c798ee6bccec9925fef485260f41baa5a`,
and clean inference executor `9d8b608ab9b4a152d67ceacc04fe49b4d7595788`.
It verifies command, working directory, task identity, timestamps, execution
logs, current runtime/input identities and the complete recorded output
inventory. It refuses pending jobs before reading native output files.

Content checks cover:

- Exactly 984,137 input genes across 78 corrected proteomes, with original
  copied filenames and bytes.
- OrthoFinder 3.1.5 native log command and successful completion.
- Exact species/sequence maps and byte-identical sequence strings in internal
  FASTAs, not just matching total sequence counts.
- Finite native graph weights and complete, unique MCL checkpoint membership.
- A rooted, node-labeled species tree containing each expected taxon once.
- All 6,006 directed native orthologue tables, including empty ones; correct
  identifiers/species ownership; and exact pair-set agreement across reverse
  orientations. Every relation must stay within one MCL checkpoint group.
- Injective accession normalization, with duplicate expanded relations counted
  separately from distinct pairs rather than hidden.

The pair-table audit holds only the two orientations of one species pair at
a time. It does not promise constant memory or a measured speedup. Top-level
per-species summary TSVs are permitted by exact filename but are not inputs
to the native pair converter or semantically checked by the pair-table audit.
The outer admission still inventories and hashes them.

## Validation

The focused 80-test suite covers scheduler/execution mutations, native table
coverage, reverse-orientation disagreement, unknown/incorrectly owned genes,
MCL boundary violations, accession collisions and input-parity helpers.
The actual pending raw job 21706 was rejected before native file access.

The content validator also ran on the retained WGD OrthoFinder outputs, using
that dataset's actual 23,870-gene/four-species scope: 5,925 MCL groups, 12
directed tables and 38,572 distinct native pairs, with zero duplicate converter
relations. This is a native-format integration check, not corrected-QfO
admission, a new biological result, or evidence that a species tree is true.

The batch wrapper requests 2 CPUs, 64 GiB and four hours on bizon and should
run after corrected OrthoFinder terminates. Native pairs, MCL clique pairs,
and the supplied-tree FastOMA workflow still need separately admitted
conversion/execution/scoring steps. Shared-host inference durations must not
be reported as dedicated matched timing.
