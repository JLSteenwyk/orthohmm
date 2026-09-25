# Recovery Batch 15 Interruption

On resumption on 25 September, fresh Slurm accounting reports native task
`22103_14` (raw job `22143`) as `TIMEOUT`, exit `0:0`, with its batch step
`CANCELLED`. It started `2026-09-24T01:14:56` and ended
`2026-09-25T12:06:17`; recorded elapsed time is `1-10:51:21` against a
one-day limit. `uptime -s` reports host boot at `2026-09-25 12:05:20`.
These records establish a terminal interrupted attempt, not 34 hours of
measured native computation or a sequence-specific BLAST failure. The exact
shutdown cause and timing have not been established.

The native log and time file are empty. The status remains `starting`, with
no native completion marker or final `hits.blast`. The partial table is
8,556,544 bytes, last modified `2026-09-24 01:17:12.650672483 -0400`.
No partial rows are admitted or reused from this attempt.

The entire original directory was retained by a non-overwriting rename from
`benchmarks/results/qfo_blast_recovery_v1/batch_14` to
`benchmarks/results/qfo_blast_recovery_v1/batch_14_interrupted_22103_14_20260925`.
Original scheduler log `benchmarks/work/qfo_blast_recovery_22103_14.log`
remains in place. File SHA256 values before preservation were:

| File | SHA256 |
|---|---|
| `hits.blast.partial` | `6e22696fd4050d15933401c1d28b4e72c9eb585083c51ea08f5795646b51dabf` |
| `preflight.json` | `65c3ef26d3c98ca45a11a0be8e9c0f0f84f5ba36f9a4055fd55a0313aa858dc9` |
| `status.json` | `65c3ef26d3c98ca45a11a0be8e9c0f0f84f5ba36f9a4055fd55a0313aa858dc9` |
| `blast.log` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `blast.time.txt` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |

## Selective Replacement

The first fourteen completed and independently admitted batches are retained.
Only interrupted index 14 (human batch 15) was resubmitted:

```bash
sbatch --parsable --array=14 benchmark_tools/results/qfo_blast_recovery_batch_20260923.sh
```

Slurm returned `22160`. Task `22160_14` is confirmed RUNNING on `bizon`,
180 CPUs, 900 GiB, one-day limit, no requeue, starting
`2026-09-25T12:14:01`. The wrapper and scientific settings are unchanged.
The execution checkout remains clean at
`1a73dc619bf50d5cc37b59f7b7820df3f4787c5d`. Its frozen preparer rechecks
the original query, full database, native runtime and command before search.
This is an explicitly selected infrastructure recovery, not an automatic
retry, best-run selection, sensitivity change or timing experiment.

## Remaining Handoff

The old admission for index 14 (`22105_14`) has an unsatisfied dependency.
Later native/admission tasks, merge `22150` and search audit `22151` remain
pending. No dependency has been bypassed. Existing admission/merge sources
pin the old array identity; they must not be used unchanged to assert that
the replacement completed under the original job ID.

Next: implement and test an explicit, narrowly scoped replacement identity
in the admission and panel/merge validation contracts, preserving original
failed-attempt evidence. Freeze updated validator code, admit the completed
replacement, then reconnect the remaining native tasks only to successful
replacement admission. Propagate updated source identities through search
admission before launching downstream BPO conversion. Do not relabel or
overwrite old successful admissions, release the unrelated held chain, or
reuse the interrupted partial table. DGX remains deferred.
