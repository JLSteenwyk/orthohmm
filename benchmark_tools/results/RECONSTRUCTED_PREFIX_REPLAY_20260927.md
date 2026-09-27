# Conda Prefix Transformations Replayed

The [replay summary](reconstructed_prefix_replay_20260927.json) strengthens
the earlier [payload audit](RECONSTRUCTED_CONDA_PAYLOADS_20260927.md): all
195 prefix-rewritten files now match payloads regenerated from their pinned
package archives, not merely the installer's recorded final digests.

For each relevant archive, verified its retained SHA256 and read its original
`info/paths.json`. The selected member's original digest, placeholder and
file mode must agree with the audited metadata. Then read and hash the raw
regular-file member and apply Conda's `replace_prefix` followed by
`replace_long_shebang`, using the actual new prefix and recorded package
subdirectory. This follows the inspected Linux path of `update_prefix`.
Predicted length and SHA256 are compared with the existing installed file.
All transformations occur in memory; no installed payload was changed.

All **195 of 195** comparisons pass: 159 binary-mode and 36 text-mode files
across 11 packages. Counts by package/mode are in the summary. The complete
member-level report and executed command are retained at
`/tmp/orthohmm-base-reconstruction-20260927/prefix_replay.json`, pinned by
the summary along with both Conda implementation files and input receipts.
All installed files and referenced inputs were rechecked after replay.

This reuses the same Conda transformation implementation as installation;
it is not an independent implementation or proof of upstream algorithmic
correctness. It does establish agreement between pinned archive inputs,
recorded installation prefix, reproducible transformation and observed bytes.
The historical weaker audit and its one Python bytecode serialization
difference remain unchanged. Generated/extra files, whole-runtime closure,
OS dependencies, rights/security and full scientific validation remain
separate requirements. Job 22337 was not modified or admitted by this check.
