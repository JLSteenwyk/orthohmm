# Separate Base Runtime Passes Integrated Fixture

The [execution receipt](reconstructed_base_fixture_20260927.json) records a
fresh base-Python reconstruction followed by the complete eight-stage
installation/inference/readback/scoring fixture. The local root is
`/tmp/orthohmm-base-reconstruction-20260927`.

## Acquisition And Installation

Acquired all 19 packages in the previously recorded base dependency graph:
51,094,707 archive bytes, each matching its retained SHA256. Five previously
downloaded exact archives were reused; the remaining 14 were downloaded from
their recorded HTTPS providers. The receipt retains every URL, filename,
version/build, original metadata identity and acquired archive digest.

An explicit file lists the local archives with their retained MD5 values;
SHA256 was independently checked before use. The existing Conda installer
was invoked offline with copies, no default packages, a dedicated cache/home,
disabled plugins and a new prefix. The first attempt rejected the reserved
destination name `base` before installation. Its logs remain intact. A
separate `python-runtime` destination then installed successfully with all
19 expected versions/builds and no additional Conda packages. A retained
implicit-default-channel deprecation warning is not an installation failure;
no global Conda configuration was changed to suppress it.

The missing Conda-owned pip dependency was handled as an explicit wheel
overlay, not silently declared satisfied in the old graph. Pip 26.2.1 was
bootstrapped by the new interpreter from the verified wheel, using isolated
mode, no index, required hashes and binary-only installation. The new prefix
passes `pip check`; its installed pip payload has 475 byte-matching files.
The base environment is deliberately a documented Conda-plus-wheel bootstrap,
not a solver-complete Conda lock or a clone of the old site-packages.

## Executed Fixture

The unchanged integrated controller was copied outside the checkout. Both
`--base-python` and `--installer-python` point to the new interpreter. It
created two fresh venvs, installed the same inference/reader locks offline,
ran native search through inferred phylogeny/reconciliation, then ran four
scientific readers and scoring. All eight stages completed once. No native
checkpoint was reused; no native fixture retry occurred.

The 16-gene fixture produced three groups and 36 ortholog pairs, four
duplications, three speciations, two marker families, one reconciled family
and two bypassed families. Root-group bytes, score fields and scientific
summary match the earlier integrated fixture after excluding its output
directory field. Synthetic fixture reference scores are integration checks,
not biological accuracy evidence; no satellite merge is exercised.

Post-run inventories match the installation reports: one bootstrap package,
11 inference packages and five reader packages. All 475 + 3,052 + 2,702 =
6,229 audited wheel payload files match their wheel bytes. This does not
audit every Conda-installed file after prefix transformation; generated
metadata/bytecode and relocated non-site wheel data exclusions remain.

The complete fixture process tree was file-traced. Neither
`/home/bizon/anaconda3` nor the original repository prefix occurs in that
trace. Acquisition and Conda installation used the original interpreter
before this traced fixture; they are not claimed independent of it. The
trace is a literal-prefix check, not an access sandbox. OS libraries remain
shared, and this is a same-host test, not cross-host restoration.

## Reproduction And Limits

The receipt binds the local `acquisition.json`, `explicit.txt`, both Conda
attempts, pip bootstrap command/lock/report, traced fixture command, copied
controller, output/readback and payload audit. Reproduction must use fresh
destinations and verify archive SHA256 values before installation. The
fixture uses the existing separately supplied native assets, reader export,
wheelhouses and test-data manifest; these are not bundled by this receipt.

No file used by full job 22337 was changed. Its dataset-scale admission is
still pending and cannot inherit success from this fixture. This work does
not establish full scientific reproducibility across hosts, controlled
timing, current security, complete runtime closure or redistribution rights.
Historical versions remain provenance rather than installation recommendations.
