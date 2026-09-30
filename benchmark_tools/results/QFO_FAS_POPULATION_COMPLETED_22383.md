# Complete Retained FAS Population Audit

The separately justified completion attempt **22383 completed with exit 0:0**
on the local Threadripper. Controller and allocation/batch accounting agree;
elapsed allocation time was 24:23, with one CPU, 16 GiB and no restart/requeue.
This is audit duration, not inference timing or a controlled tool comparison.
Null job accounting leaves CPU use and peak memory unavailable.

The original [22382 timeout](QFO_FAS_POPULATION_TIMEOUT_22382.md) remains
retained. Following the [prospective amendment](QFO_FAS_COMPLETION_AMENDMENT_20260930.md),
six immutable rows were reused and only FastOMA/OrthoMCL were freshly scanned.
No native inference, annotation, original six database joins or scientific
score/default was repeated or changed. The complete lookup was necessarily
rebuilt once because its original in-memory index was not saved.

## Result And Interpretation

The [complete eight-method table](qfo_fas_population_completed_22383/scores.md)
and [producer report](qfo_fas_population_completed_22383/report.json) retain
all fixed methods. Every eligible count matches the corrected v7 comparison;
all precomputed/missing/unannotated counts match native logs and all saved
sample classifications/lookup values pass the producer's checks. The lookup
summary matches the first pass: 59,962,787 entries, zero invalid-value entries,
46,325,666 canonical pairs and 13,637,121 canonical overwrites.

Pairs absent from the valid precomputed lookup constitute 23.58% of the
high-sensitivity OrthoHMM eligible set, 5.13% of satellite_v2, 16.86% of full
OrthoFinder and 82.05% of its sequence-only checkpoint. Across all eight
methods, this fraction ranges from 0.33% (ProteinOrtho) to 82.05%.
This describes lookup coverage, not a false-negative rate or evidence of bias.
The table's "Uncomputed" column means absent from that lookup: some such pairs
were newly scored in the small native samples. Those new scores are not used
to tighten these conservative population bounds.

For precomputed-score sum S, precomputed count P, lookup-missing count M and
eligible count N=P+M, the displayed precomputed-only mean is S/P and the
conditional full-mean outer bounds are [S/N, (S+M)/N]. They assume hypothetical
lookup-missing scores in [0,1]. They are neither confidence intervals nor
replacement native FAS scores; underlying feature-score correctness and the
sampling law are not validated by these aggregate calculations.

Satellite_v2's bounds are 0.733966-0.785282 and full OrthoFinder's are
0.605946-0.774511. These ranges overlap; the high-sensitivity/full-OrthoFinder
ranges also overlap. The previously higher native OrthoHMM sample means do
not establish a full-eligible-set advantage from these bounds. Different
prediction sets and coverage further prevent interpreting architecture
similarity alone as overall orthology accuracy. No comparator ranking follows.

## Independent Readback And Provenance

The [independent readback](qfo_fas_population_completion_readback_22383.json)
uses a separate standard-library implementation, without importing the
producer or renderer. It verifies all eight identities, conservation of
counts, exact equality of the six reused row payloads with retained originals,
native logged populations, database-command paths, seven submitted Git source
bindings and 337 current evidence identities. Rational arithmetic reproduces
the stored means and bounds within 1e-15. It does not repeat the large joins
or independently rescore the lookup. All 22 focused readback tests pass,
including rejection of partial panels, changed rows/arithmetic/provenance,
duplicate controller fields, unsuccessful jobs and changed allocations.

The [terminal receipt](qfo_fas_population_terminal_22383.json) binds submission,
producer output, readback, table and raw scheduler/accounting evidence. Thirteen
small output/controller/log artifacts were copied byte-for-byte; large databases
and lookup data remain in place. Raw whitespace is preserved. The terminal
recorder retained 54 polls without an observation error, not a resource monitor.

Current input/source/database checks before and after the completion pass do
not retroactively establish the original parser-file identity or repair its
missing final stability pass. New database hashes are retained identities,
not independent historical prediction/conversion provenance. Full eligible
exposure to complexity exclusions and dependence-aware FAS comparison
uncertainty remain unresolved. The complete publication goal is still active;
controlled timing, original TreeFam sources and final release work remain open.
