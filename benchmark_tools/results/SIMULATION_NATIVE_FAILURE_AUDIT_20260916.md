# Simulation Generation and Native Failure Audit

## Generation

Slurm array 20920 has 40 terminal tasks, all `COMPLETED` with exit `0:0`.
All 70 planned dataset truth files exist; no derived condition was
inapplicable. The pinned history verifier checked completion evidence,
generation commands, manifest provenance and output checksums for all native
runs. All 20 prespecified T/G history comparisons matched, with 437 biological
files on each side of every comparison. Copied parameter files were excluded.
The complete inventories are in `simulation_histories_20260916.json`.

These checks establish the intended matched histories and artifact integrity,
not realistic coverage of all evolutionary processes or method accuracy.

## OrthoFinder Failure Despite Exit Zero

The baseline_20261001 method task finished all three child processes with
exit zero. Both OrthoHMM configurations passed native completion and source/
input/output provenance checks. OrthoFinder did not: its native `Log.txt`
has the correct 3.1.5 version and frozen command but no run-completed marker.
The console log reports zero usable species-tree orthogroups, empty MSAs and
an error. Its graph contains `nan` and `inf` edge weights.

The eight input proteomes have 99, 97, 97, 97, 98, 98, 98 and 100 sequences,
respectively. Every sequence is 300 residues, as specified in the frozen
protocol. DIAMOND produced actual self hits and cross-species hits; the
warning about searches against a species itself does not establish that
DIAMOND failed. The native reader excludes identical-query/target self hits
by default before that warning is generated.

Calling the installed, unmodified OrthoFinder normalization implementation on
the saved species 0 versus 1 search results reproduces the numerical failure:

- 98 raw hits, with length product 90,000 for every hit.
- Fitted slope 217.07968906052554 and intercept -1072.8018499530745.
- Overflow and invalid-multiply warnings during native normalization.
- All 98 stored normalized values are nonfinite.

`diagnose_orthofinder_constant_lengths.py` reproduces this calculation using
the OrthoFinder interpreter. Its machine-readable report records the source,
input and diagnostic hashes, versions, command, warnings and counts in
`orthofinder_constant_length_diagnostic_20260916.json`. It does not patch
OrthoFinder, change sequences, generate predictions or calculate accuracy.

With identical length products, the two fitted coefficients cannot be
identified separately. The observed coefficients cause numerical overflow
in the native factorized normalization. This is a reproduced constant-length
failure case, not evidence of general OrthoHMM superiority.

## Consequences

Keep the frozen panel and its native failures as evidence. Do not silently
replace outputs, treat zero exit codes as native success, or impute failed
accuracy as zero. `validate_simulation_outputs.py` now checks native
completion/provenance and rejects nonfinite OrthoFinder graph weights.
The MCL checkpoint must not be presented as a valid sequence-only accuracy
result merely because a file exists after invalid numerical preprocessing.

A separately frozen, biologically broader panel with heterogeneous family
lengths is required for a useful general-purpose comparison. Prespecify it
before scoring or tuning on new outcomes, retain the fixed-length panel as a
stress test, and distinguish unmodified competitor results from any optional
diagnostic numerical repair. No competitor code or frozen scientific settings
have been changed. Native validation and scoring of the remaining tasks,
the additional simulation protocol, and the wider publication work remain.
