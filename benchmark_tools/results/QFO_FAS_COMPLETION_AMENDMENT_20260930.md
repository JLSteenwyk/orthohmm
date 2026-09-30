# FAS Population Completion Amendment

Prospective continuation of the fixed eight-method audit, frozen before any
FastOMA/OrthoMCL population results are inspected. The original attempt 22382
remains TIMEOUT, and its [failure evidence](QFO_FAS_POPULATION_TIMEOUT_22382.md)
is not altered or replaced. This is a separately justified completion workflow,
not an automatic retry of inference or a retry to obtain a preferred score.

## Reason And Unchanged Scope

The original 90-minute, one-CPU budget ended after six consistent partial
recounts. No count/sample discrepancy caused that failure. Repeating those
large DISTINCT joins, especially the 222,204,921-row OrthoFinder sequence-only
query, is unnecessary if their source rows and database bytes remain unchanged.
The user's goal explicitly calls for preserving scientifically reusable work.

Keep all eight corrected-v7 identities, exact native accession query, alias
rules, annotation precedence, canonical lookup overwrite order, native logged
counts, saved-sample checks and hypothetical [0,1] completion bounds. There
is no changed cutoff, reference, native FAS score, scientific default, endpoint
or selected passing-method table. All eight rows are required in the final
report, regardless of their values.

## Reuse And Evidence Gates

The fixed prior receipt SHA256 is
`7734520254562f63dbbf4fed81a2f750fd43530e5d03cd22e85888b467989461`.
Reuse exactly the six rows in its original manifest prefix, not a selectable
subset. Validate their source JSON hashes, schemas, arithmetic, database bytes,
native command database paths and original source/protocol/script Git objects.
The corrected manifest, native sample-attrition data and annotation review must
also match the original producer revision `0591bcbe`. Fail on changed bytes,
sidecars, incomplete identities, a changed query scope or a discrepancy.

Rebuild the complete lookup once with the unchanged original helper and compare
its summary to the first pass. This repeat is necessary because the original
process's compact lookup index was not persisted and is needed for the two
remaining joins and fresh sample-value checks. Do not repeat annotations,
orthology inference or any of the six completed database joins.

Freshly scan only FastOMA and OrthoMCL, using the original read-only query and
population summarizer. Recheck all eight saved sample classifications/values
against the rebuilt lookup, and all database/source/input identities before
and after completion. All eight rows, their counts, original logging, sample
checks and hypothetical-bound arithmetic must pass before `report.json` exists.

The new report explicitly says six rows are reused and two are fresh. It
preserves the original timeout, missing original final stability pass and
unestablished historical parser-file identity. Current stability must not be
described as retroactive original parser or historical prediction provenance.
The renderer refuses a selected panel, hidden reuse or relabeling this
continuation as a fresh eight-database pass. Bounds remain conditional
hypothetical ranges, not confidence intervals or replacement native FAS scores.

## Execution

One completion allocation on the local Threadripper: one CPU, 16 GiB and a
180-minute maximum; no GPU, exclusive-host claim, requeue, extension or automatic
retry. The longer allowance addresses the observed resource timeout rather
than changing scientific criteria. Stop and preserve any failure/discrepancy.
Use the existing interpreter with NumPy 2.2.6/ijson 3.5.0 and numerical-library
threads set to one; do not install/upgrade shared packages or interrupt others.
Controlled timing stays deferred and no DGX is involved.

```sh
/home/bizon/anaconda3/bin/python -B -m benchmark_tools.complete_fas_population \
  --output benchmarks/work/qfo_fas_population_completion_20260930
```

Commit and push helper, tests, amendment and allocation script before submitting
this single completion attempt. Preserve submission intent/response and the
authoritative terminal state. On success, independently read back all eight
rows and reuse provenance before rendering the complete table. Passing unit
tests or a running allocation do not admit a scientific result or publication
readiness. The original failure and partial artifacts remain retained.
