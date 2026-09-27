# Three-Way Native Prediction Comparison

Added the comparison layer for the prespecified historical, fresh retained,
and canonical QfO contrasts. It validates complete root partitions, streams
native pair/confidence tables, attributes unique pairs and changed annotations
to each side's source families, and compares explicitly rooted species-tree
clades separately from serialized bytes. Group/family relabeling is not
treated as a changed native pair prediction. Branch lengths and support are
not part of the topology equality endpoint.

The [fixture result](qfo_three_way_fixture_20260927.json) compares the three
previously checked 16-gene historical/fresh/cache-reused fixtures. All three
contrasts have identical three-group partitions, all 36 pairs/confidences,
and species-tree topology/bytes, with no affected root or pair families.
Thirty-one focused comparison/readback tests pass, including altered topology,
confidence-only differences, family attribution, relabeling and invalid pairs.

This helper is not a runtime or scientific admission gate. Full QfO outputs
must pass their scheduler/provenance checks and four readers before invoking
it. The full canonical run is still unsubmitted. Native job 22329 remains
RUNNING and readback 22330 waits on successful completion. No full-data result,
accuracy effect or benchmark-score replacement follows from this fixture.
The existing queued source inventory was not modified; this is a new helper.
