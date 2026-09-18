# Native OrthoMCL BPO And Index Checks

## Retained Probe

`probe_orthomcl_bpo_parity.py` compares the existing streaming converter with
the installed OrthoMCL 1.4 `blast_parse` function through BioPerl. The probe
does not change either converter. Four synthetic input proteins and eight m8
HSP rows exercise self hits, multiple HSPs, gaps, rounded percent identity,
zero E-values, abbreviated `e-20`/`e-200` notation, inclusion at `1e-5`, and
exclusion above the cutoff.

Both implementations produced exactly the same 225-byte BPO: six directed
pair records across three query blocks. Native `constructIDX_for_bpofile`
and `constructSE_for_bpofile` generated indexes that matched independent
Python byte offsets and inclusive one-based query ranges, including the EOF
sentinel (seven offset entries).

The saved report is `orthomcl_bpo_native_parity_20260918.json`; raw fixture,
native BPO/indexes, streaming BPO and logs remain in
`benchmarks/work/orthomcl_bpo_parity_probe_v2_20260918`. The earlier probe is
retained separately. The report records commands, tool/helper/input hashes,
loaded Perl modules, native versions and generated artifacts. Loaded modules
are inventoried after execution; this is not a hermetic production runtime
freeze or before/after proof for every loaded dependency.

## Full-File Index Validator

`validate_orthomcl_bpo_indexes.pl` independently retrieves the native Storable
indexes, streams the BPO, and checks every byte offset and query range. It
rejects malformed/nonsequential similarity IDs, missing/extra offsets,
incorrect EOF offsets, missing/extra ranges and repeated query blocks. It
prints a small JSON summary only after all checks succeed.

The BPO is streamed, but native Storable offset/range structures are loaded
in memory. This avoids an additional full JSON index export; it is not a
constant-memory index algorithm. The index validator does not validate every
BPO alignment/score field against the source BLAST table.

The installed legacy Perl initially serialized the record count as a string
after interpolation. Explicit numeric conversion fixed this schema issue;
the retained v2 probe and native regression test check the integer count.

## Verification And Limits

All 24 focused converter/index tests passed, including the installed
OrthoMCL/BioPerl parity probe. The production index script passes syntax
checking under the installed legacy Perl. Corruption tests exercise offsets,
EOF/cardinality, query ranges, native index types, record IDs and noncontiguous
queries. These tests verify the validator's rejection behavior, not production
data that do not exist yet.

Small-fixture byte parity does not prove equivalence on all m8 inputs. Before
corrected production conversion, pin the Perl/BioPerl runtime and configured
native modules, validate complete BPO content against the admitted BLAST
table, and integrate both native indexes with terminal execution provenance.
No corrected BPO job is submitted by this probe, no historical result is
overwritten, and no accuracy or publication claim follows.
