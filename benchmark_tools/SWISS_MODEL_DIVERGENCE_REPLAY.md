# Retained Model-Distance Arithmetic Replay

This prospective component reproduces distance/table arithmetic for the new
SwissTrees model-distance analysis from relocated direct inputs. It does not
repeat IQ-TREE, admission, initial scoring or bootstrap draws. The original
feature and score readers remain immutable; the replay invokes only their
frozen edge-split, rational-statistic and table-check functions.

## Build

Commit and push the tested new replay/guide/requirements before a selected build.
Build from that exact commit, not a dirty checkout:

```sh
python -B -m benchmark_tools.reproduce_swiss_model_divergence build \
  --repo . --revision COMMITTED_REPLAY_REVISION --output /fresh/component
```

The component includes 12 hash-pinned retained inputs/readers and three committed
support files, including `replay.py`. Retain the returned `REPLAY_INDEX.json`
SHA-256 outside the component before transporting it. Verify every regular file
against this anchor before executing copied code. Never overwrite an existing
component or replay output.

## Replay

Use Python 3.10 or newer with the exact requirements in the component;
Biopython 1.87 supplies the same explicitly shared parser. The existing tested
3.10.13 scientific environment is sufficient; no installation is necessary on
the current machine. The component itself does not archive wheels or an OS
runtime, and installation/availability on other hosts is not established.

```sh
python -I -B /relocated/component/replay.py replay /relocated/component \
  --index-sha256 EXTERNALLY_RETAINED_INDEX_SHA256 --output /fresh/replay
```

The command validates the complete component, then checks every one of the 183
regular evidence members against the original feature-readback identities
before extracting. It recomputes all 10,765 edge-path distances, all family
descriptors and the original median/tie bins. It replays exact-rational
macro-family precision/recall and harmonic F1 over all 54 count records, checking
every one of nine score rows, six conditional differences, both TSVs and all
15 human-table rows. The prior-adjusted counts and original outcomes remain
unchanged. The original failed R1 timing is not promoted to an eligible result.

Paths embedded in historical receipts are provenance, not fallback reads.
The replay does not call the original full `verify()` functions, which require
historical Git, scheduler, installed inference binary and logical input paths.
It needs none of those facilities for this bounded arithmetic stage.

## Limits

This is same-function replay with shared Biopython parsing, not a second
independent algorithm, biological validation, model-adequacy test, causal
mechanism, confidence interval or independent confirmation. Raw count/tree
admission is inherited. Selected inference, native OrthoHMM, alignments, species
history and failure causes are not reproduced. No timing repair, isolated
efficiency, transitive study/runtime closure, redistribution clearance,
archival DOI or publication readiness follows. Actual execution and restoration
receipts are required before claiming the replay completed.
