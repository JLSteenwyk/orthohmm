# Historical Fragment Annotation Acquisition

Exploratory, development-exposed follow-up. Aggregate and prior sequence-strata
outcomes have already been inspected. No method tuning or confirmatory claim.
A source-feasibility pilot inspected F7BEF0 before this protocol: its corrected
sequence version 3 first appears in December 2020, whereas release 2020_04
contains version 2. This pilot was selected to test corrected-input provenance,
not from a prediction error. No fragment-stratified outcomes have been joined.

## Frozen Selection Rule

Use all 563 accessions and 18 families from corrected sequence inventory
SHA-256 `912363276fa5456f11aeff9f398fce71aec8d377c8c8eda7b144105945e554f1`.
Re-extract exact sequences from its checked FASTAs, accession, taxon ID and
FASTA SV field. Never derive membership from method predictions.

For each accession retain the UniSave history response and select the entry
with matching sequence version present on 12 August 2020 (release 2020_04).
If that sequence version first appears later, select its earliest historical
entry and explicitly label it later-sequence-version annotation. Never replace
this with the current entry. Require primary accession, taxonomy ID, sequence
version, entry version and exact sequence SHA-256 to match. Failure is missing
annotation, not evidence of completeness. Do not scan alternative versions
until one produces a preferred fragment label or silently normalize residues.
Record annotation date, reviewed/unreviewed status, raw response identities
and acquisition URLs. Preserve HTTP/parser/matching failures; no overwrite.

Read structured Swiss-Prot flat files through Bio.SwissProt. Retain DE Flags
Fragment/Fragments and NON_TER/NON_CONS features separately. A protein is
annotation-positive if either signal is present. A matched entry lacking both
is unflagged, not proven complete. Taxonomy or sequence mismatch is missing,
even if the accession is shared. These are externally maintained annotation
signals, not independent experimental validation of sequence completeness.

## Outcome Display

Freeze complete acquisition/membership results before joining predictions.
Family bins: any annotation-positive member; no positive and all members
matched/unflagged; no positive and at least one missing. Preserve all families,
report dates, later-version coverage and reviewed status, and retain empty
bins. Add a baseline-release-only sensitivity display, marking later entries
missing rather than substituting their older mismatched sequences.

Show all admitted corrected methods and missing OrthoMCL explicitly. Use the
existing SwissTrees raw/2+1 per-family precision/recall, macro P/R and harmonic
F1, not pooled pairs or mean family F1. Displays are descriptive only: no new
bootstrap, p-values, subgroup significance, causal explanation or superiority
claim. Sparse/empty bins are informative limitations, not reasons to change
cutoffs or omit families. This does not resolve duplication-history evidence
or all publication uncertainty requirements.

## Sources

[UniProt entry history](https://www.uniprot.org/help/entry_history) documents
separate sequence and annotation versions and the historical archive.
[UniProt non-adjacent residues](https://www.uniprot.org/help/non_cons) describes
the incomplete-sequence feature and associated fragment flag.
[UniProt API access](https://www.uniprot.org/help/api) documents programmatic
retrieval. UniSave flat-file headers carry CC BY 4.0 attribution; retain them
with any redistributed entries. Acquisition date does not imply annotation
date or membership in the original QfO release.
