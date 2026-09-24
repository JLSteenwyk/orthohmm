# Recovery Admission 22155 Failure

Native recovery 22154 completed successfully in 11:38, including its optimizer,
refinement and repeat refinement. Independent admission 22155 failed after 1:47
when its separate refinement child received SIGSEGV. Candidate job 22156 is
PENDING(DependencyNeverSatisfied) on 22155; its dependency remains unchanged.

## Observed Evidence

The Python fault handler reports garbage collection in `read_partition`, called
at line 50 of the frozen recovery runner after `write_clusters`. The failed
child wrote `refinement_repeat.txt` but did not write its JSON completion report.
The final file is 23,875,927 bytes with SHA-256
`f4c6f1973bc9636828baf1e6d9be3f416a18fc8ada5302fb502c082495fc1811`.
Both successful native refinement outputs have identical size and digest.

A separate system-Python, standard-library-only streaming reader verified all
984,137 saved gene names occur exactly once in 390,845 groups, with no unknown
members. Input identities were checked before and after the readback. This
reader imports no NumPy, SciPy, Numba, Biopython, igraph or Leiden extensions.
[Machine-readable evidence](qfo_failed_refinement_readback_22155.json).

The reader's ten tests cover complete membership, duplicate genes within and
across groups, unknown and missing genes, wrong group count, duplicate/empty
universe, exact hashing, and CLI operation with site packages disabled.
`coredumpctl` is unavailable on this machine; the retained Python fault-handler
stack is not a native backtrace. No root cause is established. In particular,
garbage collection at the crash location does not prove a garbage-collector bug,
and byte equality does not establish safe native execution.

## Bounded Next Diagnostic

Before any revised admission, prepare one separately recorded refinement-only
diagnostic using the exact frozen native runner and saved profile/graph/numeric
inputs. Do not invoke optimization, change thresholds, disable garbage
collection, or bypass the missing child report. Use a fresh output directory,
one CPU, 64 GiB, one-hour limit, no requeue; retain inherited scientific thread
limits and enable `PYTHONMALLOC=debug` and `PYTHONFAULTHANDLER=1` for allocator
diagnostics. Record the allocation, environment difference, source/input hashes,
child return code and logs whether it succeeds or fails. Do not retry it
automatically or select among repeated successful partitions.

If the diagnostic succeeds, compare its complete membership, numeric audit and
source metadata to the existing successful refinements. If it fails, retain the
stack/output and investigate the first implicated native boundary. Either
outcome is diagnostic only, not a new seed admission or accuracy result. Any
revised admission must have a separately reviewed contract, retain failed job
22155, and require successful execution of all validation steps. No existing
failure record or blocked dependency may be relabeled as success.

This plan does not claim that allocator debugging repairs the failure. The
diagnostic has not yet been implemented or submitted. DGX remains deferred.
