# Installed Native Fixture Follow-Up

All ten tests skipped by the [full unit run](PUBLICATION_TEST_REFRESH_20260927.md)
were enabled with `ORTHOHMM_LEGACY_BLAST_SMOKE=1` and passed: ten passed,
128 unrelated tests deselected, zero skipped/failures/errors, exit zero in
21.39 seconds. The JUnit class/name set exactly matches the ten earlier skips.
Both runs remain separate records; the original skipped result is not rewritten.

The [machine-readable receipt](publication_native_unit_refresh_20260927.json)
records the command, environment override, test-source identities, current
revision, native fixture reports and JUnit checksum. Tests cover installed
legacy query/database normalization, short/regular BLAST smoke cases,
BioPerl/BPO conversion and index parity, Perl lookup isolation, checkpoint
preparation and independent rereading, and staged native OrthoMCL inference.
The staged fixture asserts twelve groups, 41 grouped proteins and 379 indexed
records. These are small fixtures, not a benchmark-dataset rerun.

Inspected fixture scope before enabling the opt-in flag. A fresh retained
temporary root, `benchmarks/work/publication_native_unit_refresh_20260927`,
keeps outputs separate from historical runs. Native tests use installed local
tools and their existing bounded subprocess checks. No production inputs,
historical outputs, scientific settings or unrelated jobs were changed.

The current source is identical to the full-unit revision across `tests/unit`,
`benchmark_tools/*.py` and `orthohmm`; scoped worktree status is clean afterward.
Together the two runs exercise all 10,739 collected unit cases on this host,
including opt-in fixtures, while preserving the original 22 frozen-source
syntax warnings. This does not establish complete cross-host integration,
independent biological accuracy, dedicated timing or publication readiness.
