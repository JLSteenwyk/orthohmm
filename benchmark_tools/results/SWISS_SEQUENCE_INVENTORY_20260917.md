# SwissTrees Sequence Descriptor Inventory

This prediction-independent inventory uses the frozen78input FASTAs and
the existing563-protein,18-family SwissTrees membership universe. The script
reads membership fields only, not method TP/FP/FN or accuracy statistics.
It verifies frozen manifest identities and each FASTA before and after
parsing. Exact accession matching is required; ambiguous matches fail.

```sh
python benchmark_tools/inventory_swiss_sequences.py \
  --counts benchmark_tools/results/qfo_swiss_comparator_counts_20260917.json \
  --prepared benchmark_tools/results/qfo_factorial_prepared_20260917.json \
  --output benchmark_tools/results/swiss_sequence_inventory_20260917.json
```

Result SHA-256:
`f4a147ff1ecdf0c8df046ab1a5c63069fef106fa627d52eaec008aae250527a5`.
All549exactly matched proteins contain canonical residues only. Their median
length is376residues and median global Shannon entropy is4.103450741bits.
Per-protein measurements, original descriptions, source filenames and
per-family summaries are retained. Entropy uses canonical-residue frequencies;
noncanonical symbols, when present, are counted separately, not silently
included as additional amino acids or assigned zero entropy.

## Missing Identity Resolution

Fourteen reference accessions have no exact input-header match:
A0A6I8Q068, B1H1F6, F6PXE7, F6QLB4, F6SX31, F6YNC0, F7AMC8,
F7BEI0, F7CLH6, F7D5A4, F7D5I9, F7DF67, Q0VGW2 and Q505H9.
They occur in APP, ASTER, HOX, NOX, POP and VATB. This is an accession
resolution gap, not evidence that those proteins are absent from the
benchmark's mapped sequence universe. The scorer may recognize aliases.
Resolve against the frozen reference mapping with unambiguous provenance
before claiming complete family-level sequence descriptors or stratifying
the full reference universe. No fuzzy mapping or descriptor imputation occurs.

## Fragment And Composition Boundaries

None of the549matched FASTA descriptions contains the literal parenthesized
Fragment/Fragments label. This cannot establish absence of actual fragments:
the retained headers are not a validated completeness annotation source.
There is no positive labeled bin for a fragment-versus-nonfragment analysis
from these exact matches. Short sequences or low entropy must not be relabeled
as fragments. Independent fragment annotations remain an unmet requirement.

Global entropy is a composition descriptor, not local low complexity,
evolutionary divergence or an error mechanism. No strata, score contrasts,
confidence intervals, method changes or new defaults were selected here.
Any subsequent composition analysis needs a frozen protocol and explicit
missingness handling before inspecting stratified scores.

Validation:12unit tests cover entropy, unknown-only sequences, invalid input,
literal fragment labels, missing matches, duplicate/shared accessions,
unexpected header format and changed source files.
