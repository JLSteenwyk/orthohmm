# Incomplete FAS Population Recount

The single prospective analysis attempt **22382 timed out** at its original
one-CPU, 16-GiB, 90-minute allocation. The
[terminal receipt](qfo_fas_population_terminal_22382.json) retains controller,
accounting, log, source bindings and all six partial row outputs. No scientific
inference was rerun, allocation extended, source changed or retry submitted.
The [protocol](QFO_FAS_POPULATION_PROTOCOL_20260930.md) remains unchanged.

## Outcome

The complete frozen lookup was streamed: 59,962,787 entries, zero invalid-value
entries, 46,325,666 valid canonical pairs and 13,637,121 canonical overwrites.
Six partial rows report matching native precomputed/missing/unannotated counts
and saved-sample lookup checks: both OrthoHMM configurations, both OrthoFinder
outputs, SonicParanoid and ProteinOrtho. There is no completed FastOMA or
OrthoMCL row and no final `report.json`. The log does not establish how far
FastOMA processing progressed before termination.

**This is not a completed eight-method recount.** No selected-method score
table or full-population admission is published. The six partial checks do
not replace the omitted final all-input/parser/source stability pass. The
original native eight-method FAS scores remain unchanged.

## Verification

The fresh resume recorder retained the terminal controller state without
observation errors; the earlier intentionally interrupted recorder remains
historical partial evidence. The controller reports `TIMEOUT`, reason
`TimeLimit`, exit `0:15`, zero restarts and elapsed `01:30:02`. Accounting
reports allocation `TIMEOUT` with exit `0:0`, and batch `CANCELLED` with exit
`0:15`. These distinct raw records are preserved; allocation exit zero must
not be mistaken for successful analysis.

Post-timeout readback verified all six retained database hashes, count
conservation, immutable comparison eligibility, native logged count equality
and stored hypothetical-bound arithmetic. It did not repeat the DISTINCT
joins or independently rescore the lookup. The three submitted source,
protocol and script pins match exact Git objects at `0591bcbe`; recorder
source pins also remain unchanged. Twelve small raw artifacts are retained
beside the receipt, not large databases or a reconstructed complete report.
Raw producer whitespace is preserved byte-for-byte, including the controller's
trailing space. The complete-panel renderer also rejected these actual six
rows when they were presented with a claimed-complete status in a negative
guard check; no table file was generated.

The current accounting collector is `(null)`: zero recorded CPU time and
blank peak memory are unavailable measurements, not zero resource use.
Neither the allocation duration nor this shared-host audit supports a
comparative timing claim. Database hashes are new retained identities, not
proof of historical prediction conversion. Conditional hypothetical completion
bounds are not native FAS scores, confidence intervals or sampling-law evidence.

## Next Work

Preserve this failed attempt and its partial evidence. Any future completion
workflow needs a separately justified prospective amendment; there is no
automatic retry or selection of passing methods. Reuse defensible completed
work where possible rather than repeating expensive joins without a reason.
Controlled timing remains deferred on the Threadripper. Full QfO uncertainty,
original TreeFam inputs and the final publication package remain unresolved.
