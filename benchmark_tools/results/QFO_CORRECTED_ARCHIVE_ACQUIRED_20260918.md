# Corrected QfO Archive Acquired

Acquisition job21687 completed0:0 on bizon in1:46:26. It resumed the retained
79,556,608-byte partial file; curl reported HTTP206 and2,569,109,590new bytes.
Completion was logged at2026-09-18T04:18:29Z. The resulting file has exactly
2,648,666,198bytes and SHA-256
`29a0f54e4af7d6bdbfe28923efd2f3633e844b0c49b0a9e13d2b08f6a2d66009`.
A separate post-completion sha256sum invocation reproduced this value.
This is a locally measured checksum, not a publisher-provided checksum.

Source:
[EBI corrected QfO 2020 archive](https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/previous_releases/qfo_release-2020_04_with_updated_UP000008143/QfO_release_2020_04_with_updated_UP000008143.tar.gz).
Destination:
`benchmarks/work/qfo_corrected_source_20260917/QfO_release_2020_04_with_updated_UP000008143.tar.gz`.
The resume request sent `If-Match: "9ddf7056-5b478e2a3e88c"`; this records
the request condition, not a claim that response headers were archived.

Acquisition script SHA-256:
`8f0c329197a71ff2443ce8a03b13fd20332d037f525bff601ef7f58fe9362790`.
Retained `acquire_21687.log` SHA-256:
`b1b5bef7f7433d6e0eadecf49d50dd61b49ba3a3b6b5c1d35700e7682fc9c56d`.
Neither the archive nor large raw data are committed.

## Canonical Comparison

Job21688 completed0:0 in1:47 using pinned executor
`a1fbffcc8cb5b4fdabca18843d79d94b0baadaef` and Biopython1.86/Python3.10.13.
Its report `qfo_corrected_archive_comparison_20260917.json` has SHA-256
`72cd351a835ee445a64388c2e292919d8d6b70b79afe48e6c2d696dde5f37ba6`.

- All78canonical proteomes are present; only `UP000008143_8364.fasta` differs
  from the original frozen archive inputs.
- All984,137numeric reference identities are uniquely covered, with zero
  unmapped accessions and zero missing identities in every reference species.
- All14previously missing SwissTrees accessions are recovered.
- The archive stream was read through gzip EOF.

This closes the acquisition/canonical-mapping portion of the conditional
protocol, not the native sequence compatibility or rerun authorization gates.
Job21689 is performing the separate native-sequence comparison. Extraction,
direct staged inventory and an immutable inference manifest remain required.
The original input directory, scores and running experiments are unchanged.

## Staging Representation Fix

The production comparison revealed a fixture mismatch: the archive auditor
emits a78-species dictionary containing zero missing counts, not an empty
dictionary. The stager incorrectly rejected any nonempty dictionary.
It now requires78nonempty string species keys and integer-zero counts.
Positive, negative, Boolean, string and null counts, and an empty dictionary,
are rejected. Both reports must still agree and the explicit missing-ID list
must remain empty. This fixes schema interpretation, not the scientific
zero-missing criterion. All59staging/archive/direct-inventory tests pass.
