# External Phylogeny Tools Relocated and Traced

The [installed inferred-phylogeny fixture](PUBLICATION_FROZEN_PHYLOGENY_INSTALL_20260926.md)
has now completed with private local copies of MAFFT and FastTree under
`/tmp/orthohmm-external-tools-20260926`, outside their original installation
directories. The [result](publication_relocated_phylogeny_tools_20260926.json)
has SHA-256 `2af559b53600e6df745ad06eb871acbb903fc9880a5627a5139351b77fd4e7ef`.

## Relocation Boundary

The MAFFT launcher contains an absolute compiled-in helper-directory default.
Copying the launcher alone is therefore insufficient to establish relocation.
The probe preserves the launcher bytes and sets `MAFFT_BINARIES` explicitly to
the copied helper directory; it does not patch the installed script. PATH is
restricted to `/usr/bin:/bin`. Explicit CLI paths select the copied MAFFT and
FastTree entrypoints.

The 42 copied files comprise two tool entrypoints, all 34 regular files in
MAFFT's helper directory, two available MAFFT notice files, and four synthetic
input FASTAs. Byte hashes and permission modes match the source copies.
Symlinks/nonregular helper files are rejected. Original installations remain
unchanged. These are private local copies, not newly uploaded binary artifacts
or a cleared redistribution bundle. FastTree's installed directory has only
its executable; upstream source/notice acquisition remains a release task.

## Observed Execution

One fresh pipeline run, with no checkpoint reuse, returned exit zero under
`strace -f -qq -s 4096 -e trace=%file`. The 3,155,753-byte trace has SHA-256
`a2568e079385ffd0a5e129026513ab1a0a8ad592cb0126cd1488032f82b83937`.
Neither original MAFFT nor original FastTree installation prefix appears in
the trace, including failed path lookups. Successful copied-helper executions
are visible for `version`, `countlen`, `replaceu`, `tbfast`, `dvtditr`, `restoreu`
and `f2cl`, as well as the copied MAFFT/FastTree entrypoints.

The inferred pipeline again covers all 16 genes in three root HOGs, infers a
species tree from two marker families, reconciles one duplicate-containing
family and produces 36 cross-species pairs. All five selected primary output
files are byte-identical to the preceding installed run:

- Final orthogroup partition.
- Root-HOG table.
- Pairwise ortholog table.
- Rooted inferred species tree.
- Reconciled gene tree.

All 270 input/source/copy/output/log/trace record entries were rechecked.
Seventeen related tests pass, including seven new byte/mode, symlink and trace
checks. The trace gate rejects even unsuccessful access to an original tool
prefix; empty or non-executable traces fail. All copied files remain local;
the committed report contains their identities, not executable bytes.

```sh
/usr/bin/python3 -S -m benchmark_tools.probe_relocated_phylogeny_tools \
  --repo . --output /tmp/orthohmm-external-tools-20260926
python -m pytest -q tests/unit/test_probe_relocated_phylogeny_tools.py \
  tests/unit/test_verify_frozen_phylogeny_install.py
```

Existing output directories are refused. A separately justified repetition
requires a new destination. The child process group has a 180-second timeout
and no automatic retry. Inputs, scientific flags and installed wheel are
unchanged; only external-tool paths, helper override, PATH and output paths
differ from the prior command.

## Limits

This is a same-host external-tool relocation, not an isolated filesystem or
complete runtime transfer. The Python/OrthoHMM installation, system shell,
utilities and shared libraries remain on the original host. Literal trace
prefix checks are not proof against every path alias or a complete dependency
closure. Other datasets can exercise additional tool paths. The tiny synthetic
fixture establishes no biological accuracy, general numerical equivalence,
hardware portability or comparative timing advantage. File-level license and
notice review still precedes any public archive of these binaries. No DGX
access, benchmark score update or scientific-default change occurred.
