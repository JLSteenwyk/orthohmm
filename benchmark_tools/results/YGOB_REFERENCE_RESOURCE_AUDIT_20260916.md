# YGOB Reference-Resource Audit

Reviewed before YGOB inference results or accuracy scores were inspected.
Scope: the frozen YGOB validation methods, not every historical competitor,
every upstream annotation, or every executable dependency.

## Reference Construction

YGOB's foundational paper describes manual homology curation using sequence
similarity and genomic context. This provides syntenic evidence beyond the
sequence-only inputs supplied to our methods, but not sequence-independent
ground truth. The paper describes the original seven-species resource;
it does not document the complete revision history of every v7 pillar.
[Byrne and Wolfe, 2005](https://pubmed.ncbi.nlm.nih.gov/16169922/).

The acquired v7 README identifies the August 2012 release and explicitly
states that the two post-WGD positions are arbitrary, not A/B tracks.
We therefore retain the frozen curated-group-recovery endpoint, not an
A/B-resolved orthology claim. Removing Saccharomyces from inference does not
remove its historical contribution to curation of other species' pillars.
Source: `benchmarks/work/independent_ygob_v7/README`, obtained from the
[official v7 archive](http://ygob.ucd.ie/data/v7-Aug2012/).

## Frozen Inference Inputs

| Method | Observed source of biological input | Reference-label access in reviewed workflow |
| --- | --- | --- |
| OrthoHMM high_sensitivity | Prepared 16-species FASTAs; profiles built from inferred sequence groups | No reference-group argument passed into inference |
| OrthoHMM satellite_v2 | Same FASTAs; species tree inferred internally from inferred single-copy families | No supplied YGOB tree, pillar labels, or synteny tracks |
| OrthoFinder 3.1.5 full | Copies of the same FASTAs; DIAMOND databases built from species FASTAs | Fresh `-f` run, not assignment against a supplied annotated reference collection |
| OrthoFinder sequence checkpoint | MCL checkpoint from that same full run | No separate label-driven inference |

Evidence inspected:

- `run_ygob_validation.slurm` passes only prepared FASTAs and fixed settings
  to the inference commands. Its preflight reads the reference file solely
  to verify its checksum, not to select genes or change inference parameters.
- At frozen OrthoHMM commit `7f3a9e4`,
  `search/profile_expansion.py:63` loads profiles' sequence database from
  the input FASTAs and requires exact agreement with the inference gene table.
  `phylogeny_pipeline.py:1061` takes the internally inferred species-tree
  branch; the manifest describes that source at line 1071.
- Installed OrthoFinder `run/run_commands.py:120` constructs search databases
  from each species FASTA. `run/config.json:59` defines the selected DIAMOND
  database/search commands. These are input-derived databases, not evidence
  of an external orthology-reference lookup. The
  [official tutorial](https://orthofinder.github.io/OrthoFinder/tutorials/beginner-tutorial/)
  likewise describes proteome FASTAs as the input to a fresh analysis.
- The frozen production worktree has no tracked differences in `orthohmm`
  or `benchmark_production.py` relative to its pinned commit at this audit.

This is a command/source audit, not a system-call trace or a proof that all
transitive dependencies are incapable of opening other files. Substitution
matrices and evolutionary models encode previously estimated biological
information; lack of a labeled database argument does not imply that every
model parameter was learned without related sequences.

## Known Dependence

`ygob_homology_screen_20260916.json` records qualifying development hits for
71,714 of 83,404 proteins (85.98%) and 6,952 of 10,250 pillars (67.82%).
These counts demonstrate family overlap, not transfer of curated labels.
Only the best qualifying hit is retained, so the per-dataset best-hit counts
are not exhaustive overlap counts. A negative screen does not establish
family independence.

The candidate-selection audit already records substantial Saccharomyces
exact-sequence overlap and the resulting whole-genus exclusion. Other
conserved sequences and shared annotation ancestry remain possible. QfO and
OrthoBench informed method development; that exposure remains a limitation
even though YGOB labels were not used for this frozen configuration.

## Interpretation Decision

The current evidence supports evaluating the prespecified novel-taxon,
curated-group transfer experiment. It does not support calling it
family-disjoint, wholly independent of sequence-based curation, or proof of
generalization to arbitrary datasets. No groups, taxa, thresholds, or endpoints
are changed after this audit. Keep the existing accuracy gate until inference
completion, native-output verification, and scoring checks also pass.

Before making stronger independence claims, obtain per-pillar curation history
and trace shared annotation/reference ancestry where available. That stronger
claim is not needed for the bounded transfer experiment. Dataset redistribution
permissions and archival terms remain unresolved and must be settled before
redistributing the raw snapshot; this audit does not grant permission.

## Inspected File Checksums

Paths below are relative to the repository except installed OrthoFinder files,
whose prefix is `SOFTWARE/orthofinder_3.1.5/lib/python3.12/site-packages/orthofinder`
on the shared volume. Source hashes preserve the reviewed implementation even
if the working repository later changes.

| File | SHA256 |
| --- | --- |
| `benchmark_tools/run_ygob_validation.slurm` | `aa38a0c700961cff280b103fb05738943bc8c3aa3fda073fdef2a94c893e9038` |
| `benchmarks/work/independent_ygob_v7/README` | `162fb37cf0f7a50b44fae47f0d2c2f9b126b728660b65d9bee8d3175af8317f6` |
| Frozen `orthohmm/search/profile_expansion.py` | `1025cfab4f80e3a85d6eaa3e768ae5e0e25dedb48fdba9deb08cb6e6aa7bb42a` |
| Frozen `orthohmm/phylogeny_pipeline.py` | `44e00316546b5df78354badd2b1a6bb595b685e98b86f686f1a32543d7f15e4f` |
| Installed OrthoFinder `run/config.json` | `10a10a93262c9676865f03ba05b2b77d48285c312dca435d2458939894a474b3` |
| Installed OrthoFinder `run/run_commands.py` | `10aeb00c3affb31175a337478957807b79b153b9b91597b29e4dccb782d8d387` |
