# Immutable Package Archive Output Guard

The archive helper formerly accepted a fresh output inside its verified input
directory. Opening that output changed the exact inventory and could include
the archive itself in the payload walk. A path through a directory symlink had
the same effect. This is a packaging defect, not a scientific scoring change.

Four new copied-CLI fixtures reproduce it under normal/optimized isolated,
site-disabled Python with direct/symlink-aliased output paths. Before the fix,
all four fail because archive creation incorrectly succeeds. Reject the
resolved output beneath the resolved immutable input before verification or
opening the output. After the fix, all thirty package cases pass in2.09s
(JUnit2.024s), including unchanged original inventory after rejection.

Retained local JUnit evidence:

| Receipt | Tests | Failures | SHA256 |
| --- | ---: | ---: | --- |
| `benchmarks/work/publication_package_output_guard_before_20261004.xml` | 4 | 4 | `6de08aa6d676f29da5b9d88a8c5f6a3b4d379508502ef80c1a4fa24595f4eb68` |
| `benchmarks/work/publication_package_output_guard_after_20261004.xml` | 30 | 0 | `ce1447867272f19c6beff88633600f894e19e6b92c6cae6592d823d74f405094` |

Fixed helper11794 bytes, SHA256
`948efb04abfa84aea1c6822d9145592a15d54416838245f0c1bf569f5d81bbab`;
test source8410 bytes, SHA256
`bcec808b2bb46004b2c983bf5ca292b17668e24626d95bdfefebc2708d70712b`.
These tests use temporary synthetic packages, not the real study archive.

The retained [rc1 execution](publication_package_execution_20261004.json)
wrote its archive outside the immutable package and is unaffected. Preserve
its original helper/index/archive bytes and historical26-case receipt. The
current guard applies to future packaging; do not rebuild unchanged rc1 or
rerun successful native inference solely for this repair.
