# Allocator Diagnostic 22158

The single planned allocator-debug scientific execution failed. Slurm reports
FAILED 1:0 after 33 seconds; its only refinement child exited -11 (SIGSEGV).
The preceding 22157 attempt had stopped in wrapper preflight and launched no
scientific child. Both failures remain retained.

The child again reports garbage collection in the frozen `read_partition`
function at the post-write readback call. `PYTHONMALLOC=debug` did not report an
allocator violation before the segmentation fault. This neither exonerates
native extensions nor proves a CPython garbage-collection defect.

The output has the same SHA-256 as both successful native refinements and the
failed admission's output:
`f4c6f1973bc9636828baf1e6d9be3f416a18fc8ada5302fb502c082495fc1811`.
An independent streaming readback verified 984,137 unique genes in 390,845
groups. All 984 recorded diagnostic input/source references were rechecked.
No child completion JSON exists, so neither byte identity nor coverage
substitutes for a successful validation execution.

[Failure audit](qfo_cpm_allocator_failure_22158.json) and
[submission receipt](qfo_cpm_allocator_submission_22158.json) preserve exact
source, scheduler, input and output identities. No optimizer was invoked.

## Consequences

Recovery admission 22155 remains failed. Candidate job 22156 remains blocked;
no high-CPM accuracy result may enter the parameter panel. Low-CPM and threshold
results remain independently admitted; missing high-CPM uncertainty remains
explicitly missing with the full planned multiplicity correction retained.

No automatic retry or allocator-mode search will be performed. The next useful
debugging evidence would be a native backtrace or a smaller reproducer, not
another attempt to obtain a successful partition. `gdb` is installed, but no
Python core report was found in `/var/crash`; the system routes dumps through
Apport and `coredumpctl` is unavailable. No unrelated crash reports were read.
No kernel segfault/OOM entry was returned for the diagnostic's time window.
These observations do not establish hardware health or a root cause.

Further native reproduction needs its own bounded diagnostic design and must
remain separate from scientific admission. The broad publication goal remains
active, including completion of OrthoMCL recovery, defensible error analyses and
publication packaging. DGX remains deferred.
