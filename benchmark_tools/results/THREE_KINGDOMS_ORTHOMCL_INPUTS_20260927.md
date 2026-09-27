# Three Kingdoms OrthoMCL Input Content

The earlier [method-input audit](THREE_KINGDOMS_METHOD_INPUTS_20260918.md)
verified OrthoMCL's recorded staged-input hashes but did not trace its native
per-species copies or merged search FASTA. The new
[machine-readable readback](three_kingdoms_orthomcl_inputs_20260927.json)
closes that retained-file content check without inference or scoring reruns.

All twelve native `orthomcl/data/*.fa` copies match staged input bytes.
Their 443,217 unique gene identifiers and sequence-length/SHA256 inventories
match the staged panel exactly. The retained `Sep_7/tmp/all.fa` has the same
443,217 identifiers and identical sequence content: no missing, added or
changed sequences. FASTA descriptions and wrapping are not part of sequence
identity. The `all.gg` genome map assigns every gene to its original species,
with no duplicate memberships, missing species or permuted species labels.

The original native log names this merged FASTA in the formatdb and BLAST
commands and reports twelve genomes and 443,217 sequences. The recovered-run
log names the same `Sep_7/tmp/all.gg`. Both logs, the original input checksum
manifest, all inspected FASTAs and the previous staged-source report are
hashed in the receipt and rechecked after parsing. These retained records are
not immutable proof of bytes consumed by historical processes. BLAST database
internals, query transformations, recovered BPO identity and full downstream
execution provenance are outside this particular audit.

Twenty-four focused content/species-mapping/source-audit tests pass, including
changed/missing/extra sequences, wrapped FASTA identity, duplicate or empty
FASTA records and invalid genome mappings. No input, old audit, score or
scientific default was changed. The older high-sensitivity input-hash gap
remains; this result does not prove historical uniform-input equivalence for
all methods or establish proteome-wide accuracy from the BUSCO endpoint.

```bash
python -m benchmark_tools.audit_three_kingdoms_orthomcl_inputs --repo . \
  --output /absolute/new/three-kingdoms-orthomcl-inputs.json
```
