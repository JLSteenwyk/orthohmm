# Corrected-Release Descriptor Update

The initial input-only preparation stopped before writing a completion
report because its requirement that all549 old descriptors remain identical
was too strong. The corrected Xenopus release changes three sequences among
those549, not only the14 previously absent SwissTrees accessions. A fourth
record has a header-only protein-existence-level update. This is consistent
with the previously retained native sequence audit, not a new score effect.

Explicitly retain these cross-release differences:

- F6PXR2:344 to300residues, sequence version2 to3.
- F7BEF0:1033 to351residues, sequence version2 to3 and description changed.
- F7C949:1467 to451residues, sequence version2 to3.
- A0A6I8Q293:description PE=4 to PE=3, unchanged sequence descriptors.

Require the other545 old records to remain identical. Verify the three
corrected sequence hashes against their native database identities in the
frozen original-input sequence audit. Keep all old/new changed fields in
the new inventory. Do not silently discard old records or modify historical
scores. This supersedes only the original protocol's all549-unchanged check;
all composition definitions, binning rules, statistical contrasts and
interpretation limits remain unchanged. No corrected stratified outcomes
have been calculated or inspected.
