# Corrected OrthoFinder Pair Conversion

`prepare_qfo_corrected_orthofinder_pairs.py` converts the independently admitted
corrected OrthoFinder run into two separately named, unscored QfO submissions.
The full method uses only native inferred ortholog relations. The sequence-only
diagnostic uses cross-species cliques from the pre-phylogenetic MCL checkpoint.
Neither uses root-HOG cliques as a substitute for full native predictions.

## Preconditions And Accounting

Require successful terminal admission on two CPUs/bizon, its report hash and
exact source from clean executor `33e2310d1095f64e77cf068ef2a6fccd4b32a8d4`.
Require corrected 984,137-gene/78-species scope, the frozen primary plan and
scoring environment, and the admitted 6,006 directed-table audit. Rehash the
admitted evidence before and after conversion. The native admission establishes
execution/input integrity; conversion does not rerun native inference.

Native conversion reads one lexically selected orientation per species pair,
checks its raw/unique/duplicate counts against admission, and emits distinct
canonical accession pairs. Species ownership plus injective accession mapping
prevents the same canonical pair from arising in different species pairs.
Duplicate native relations are explicitly reported, not treated as independent
observations. Native conversion preserves the existing converter's pair set.

MCL conversion restores the complete unique input partition, groups genes by
species within each cluster, and emits cross-species products only. An
independent species-size formula checks its total:
`sum((group_size**2 - sum(species_size**2)) // 2)`.
Singleton and single-species clusters contribute no pairs but remain subject
to complete input-coverage validation.

Both methods require zero reference-mapping loss and unchanged evidence before
renaming partial outputs. Exceptions retain a failed result record and partial
files. Existing destinations are never overwritten. Distinct participants are
`qfo_corrected_orthofinder_full` and `qfo_corrected_orthofinder_sequence_only`.
Output status is `corrected_orthofinder_pairs_prepared_unscored`; neither score
admission nor publication readiness is implied.

## Verification

The 62-test focused suite includes conversion semantics, native duplicate
removal, exact agreement with existing native/group converters, malformed
partitions, aliases, changed audit counts, pending admissions, successful
report lifecycle, no-overwrite behavior, mapping-loss failure, and mutation
during conversion. Lifecycle fixtures mock upstream admission/runtime identity;
pair writing, checksums and mapping filters execute normally.

Retained four-species WGD outputs also passed both actual conversion functions:
38,572 distinct native pairs (exact pair-set agreement with the existing native
converter) and 54,288 distinct MCL clique pairs, with no duplicate outputs.
This is a conversion integration check, not a new accuracy analysis or proof
that corrected-QfO outputs exist.

The batch wrapper runs the two conversions sequentially, two CPUs/64 GiB per
task, with a twelve-hour ceiling on bizon. It must depend on successful native
admission job 21731. Corrected six-endpoint scoring needs its own subsequently
frozen adapter and independent score admission; original-release scores must
not be transferred to these outputs.

## Frozen Queue

Converter executor `benchmarks/work/publication_qfo_corrected_orthofinder_pairs_v1`
is detached at `aa8da7800c4801684726151ad249da7c82a3b88d`. Array 21733
contains task 0 (full native) and task 1 (MCL diagnostic), with concurrency
one. Scheduler inspection confirms the `afterok:21731` dependency, two CPUs,
64 GiB and twelve-hour ceiling. Both tasks remain pending; no corrected pair
output is claimed by this queue record.
