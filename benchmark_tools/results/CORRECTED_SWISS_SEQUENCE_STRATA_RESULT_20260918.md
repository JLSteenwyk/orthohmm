# Corrected SwissTrees Sequence Inventory

## Complete Coverage

The [input-only protocol](CORRECTED_SWISS_SEQUENCE_STRATA_PROTOCOL_20260918.md)
was committed/pushed as a52feb8 before extraction. The first attempt stopped
because four of549 old descriptors changed with the corrected Xenopus release.
The [explicit cross-release update](CORRECTED_SWISS_DESCRIPTOR_UPDATE_20260918.md)
was committed/pushed as f26c9e8 before the successful extraction. No threshold,
binning rule or statistical contrast was changed, and no accuracy outcome
was joined to these features.

All563 SwissTrees reference proteins now have exact accession matches,
including the14 absent from the old input inventory. All545 unchanged records
match exactly; three changed sequences match the native reference database's
sequence hashes, and the remaining change is the recorded PE header update.
All78 corrected FASTAs were checked before/after feature extraction.

The [machine-readable inventory](corrected_swiss_sequence_strata_20260918.json)
is351883bytes, SHA-256
`912363276fa5456f11aeff9f398fce71aec8d377c8c8eda7b144105945e554f1`.
It retains per-protein descriptions, lengths and entropy, family memberships,
all changed fields, source identities and the exact frozen strata.

## Frozen Descriptive Strata

The median family normalized-entropy cutoff is0.9540164297743712. All
proteins meet the eligibility criteria; there is no missing-entropy family.

| Primary bin | Families |
|---|---|
| Lower entropy,9 | ASTER,BAR,CASP,CITE,HOX,MAPT,PSEN,RPS,SERC |
| Higher entropy,9 | APP,BAMBI,Clusterin,GH14,NOX,POP,SUMF,TRFE,VATB |

No protein satisfies the frozen entropy<0.8 concentrated-composition flag;
do not redefine that threshold to populate a bin. The lower primary bin
therefore means relative composition differences, not established low complexity.
Seven families(APP,BAR,CITE,HOX,MAPT,PSEN,VATB) have a sequence shorter than
half their median length; eleven do not. No sequence has noncanonical residues
or an explicit Fragment/Fragments label in the retained description. None of
these observations establishes biological completeness or truncation.

The planned corrected-release primary contrasts/27-endpoint conditional
bootstrap remain unexecuted until the required corrected predictions and
counts are admitted. Historical scores are not reused as corrected outcomes.
Independent fragmentation evidence, evolutionary divergence annotation and
other endpoint uncertainties remain separate requirements.
