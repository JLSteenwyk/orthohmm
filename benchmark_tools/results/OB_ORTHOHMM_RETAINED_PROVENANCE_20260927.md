# Historical OrthoHMM OrthoBench Provenance

This [machine-readable audit](ob_orthohmm_retained_provenance_20260927.json)
supplements the [partial eight-method register](OB_PROVENANCE_REGISTER_20260926.md)
with the actual historical OrthoHMM run records. It does not borrow input,
runtime or memory evidence from the fresh installed reproduction.

## Prediction And Source Binding

The high-sensitivity score binds exactly to the `strict_profiles_refined`
output in the August 28 replay receipt, not another intermediate partition.
The cached hit file and all four recorded replay outputs match their recorded
bytes and SHA256 hashes. The replay has an input-directory command but no
per-FASTA historical checksum manifest; that input-history gap remains.

The September 2 phylogenetic score binds exactly to the root-HOG output in
its native harness manifest. All twelve recorded FASTAs and eight outputs
(including edges, clusters, root groups, native pairs, reconciliation summary,
phylogeny provenance and log) match their recorded bytes and hashes. The native
metrics report successful completion, and epoch subtraction agrees with wall
time to the recorded precision. These are retained-record checks, not immutable
attestation of historical process consumption.

All 33 recorded source entries match the corresponding Git blobs at their
recorded commits: the replay script at `37c89fdd39e1b28aeccdfcfc46c1ba80b47209f0`
and 32 phylogeny-run source entries at `bb27197e4a2f5a21139af5d6f39d8db42e483fc3`.
Both historical worktrees were recorded as dirty; matching the explicitly
listed sources does not prove that unlisted dependencies or external tool
payloads matched those commits. No historical source is silently replaced
by its current worktree version.

## Resource Scopes

| Historical row | Recorded wall seconds | Memory evidence | Scope |
| --- | ---: | --- | --- |
| High sensitivity | 319.467197 | 5.560291 GiB peak process RSS | Cached-hit downstream replay; initial search excluded |
| Phylogeny satellite_v2 | 3274.102675 | 12,733,403,136 bytes sampled sum of process-tree RSS | Full inference on shared host; scoring excluded |

The phylogenetic metrics record 73,322.972274 user CPU seconds and
3,001.015671 system CPU seconds. Stage records and exact commands are retained
in the audit. Shared pages may be counted repeatedly in summed RSS; it is not
unique physical memory and differs from the replay's process-RSS measurement.
These durations have different scopes and must not be divided into a speedup
or used as matched-resource evidence. The full-inference wall time is not the
9,438-second fresh installed reproduction and is not replaced by it.

Nine focused manifest/timing tests pass. The real audit checks 29 file records
before and after extraction, plus 33 historical Git source payloads. No scores,
defaults or historical files changed. Full transitive environment provenance,
the replay FASTA history and controlled comparative resources remain incomplete.

```bash
python -m benchmark_tools.audit_ob_orthohmm_retained --repo . \
  --output /absolute/new/orthohmm-ob-provenance.json
```
