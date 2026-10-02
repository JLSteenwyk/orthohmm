# Public Phylogeny Artifacts without Installed-Tool Prerequisites

The copied-source workflow now acquires the exact historical MAFFT/FastTree
source, binary and notice artifacts without supplied tool installations.
It downloads seven files/2,753,786 bytes, verifies their frozen identities,
and leaves every file non-executable. Nothing is compiled or run.
This is a step toward executable reproduction, not an installed native
runtime, new biological result, controlled timing panel or publication release.

## Source and Tests

Source checkpoint `c94f3a96f096a5c8ae7b388b6373dc2ffb1f6c9c` adds
`acquire_publication_phylogeny_tools.py` and focused tests. The shared
bounded downloader gains an explicit keyword-only host policy, preserving
its default base/wheel provider set and redirect handler. The tool caller
has its own restricted redirect handler; it does not broaden other workflows.

All **259 tests pass in 26.12s**, zero failures/errors/skips, in the four-module
tool/base/wheel/source panel. Retained JUnit is 38,977 bytes, SHA-256
`dcf87a894d79a33c1ab982ecd103bbb918053de2c1b1e717615c5cd528a7ad15`.
The JUnit suite records 26.073 seconds; CLI wall time is diagnostic only.
Tests exercise changed inputs/payloads, malformed URLs/redirects, response
size/hash/encoding bounds, fresh canonical output/acknowledgement, executable
permissions, later-download failure preservation and refusal to retry.
No new scientific default, historical lock, dataset, endpoint or score changes.

## Frozen Artifact Inventory

| Artifact | Bytes | SHA-256 |
| --- | ---: | --- |
| MAFFT source archive | 758,305 | `2876f4adc1a2de4ed206bc40896763bf208bf1a02bda52f8bfdd91cf52d73e4a` |
| FastTree binary | 1,496,928 | `55a9d997813aae2208bd4c2081bfa690e0ecdba2d6c491805d8689415c43e38e` |
| `FastTree.c` | 395,674 | `975202a6b74c9996af871404ff043bb2152edcbda539035662514bc12d1f3431` |
| FastTree `LICENSE` | 35,149 | `3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986` |
| FastTree `README.md` | 662 | `7575f311e6098306079988c9a4933ce21c8c7b8266f655d190b49d85dc539a7b` |
| FastTree `ChangeLog.txt` | 17,271 | `cb9388cc08330a90417571e708cdb0e2418bb59a67863bdf7772196e5b19ffe7` |
| FastTree `index.html` | 49,797 | `6b8e5747a9e959127fde76e52636e35dc09105e6817893606f2ee7f55a856326` |

The [official MAFFT source page](https://mafft.cbrc.jp/alignment/software/source.html)
lists the pinned 7.525 with-extensions archive. Its exact size/hash match the
[previous MAFFT build evidence](publication_mafft_build_20260926.json).
FastTree files come from [the immutable upstream revision](https://github.com/morgannprice/fasttree/tree/29c5e62fbcd93230ee325f9c6a17b81f00e3c72a)
and match the [previous FastTree acquisition](publication_fasttree_acquisition_20260926.json).
Those earlier receipts are read only for post-acquisition comparison, not
required by the new controller. No old installation or original artifact is
an acquisition input.

Only HTTPS, credential-free/query-free URLs on `mafft.cbrc.jp` and
`raw.githubusercontent.com` are accepted, including redirects.
Downloads are bounded by exact expected sizes and SHA-256; response encoding
and declared lengths are checked independently. Partial/failure artifacts
remain with no retry. There is no automatic upstream version selection.
The socket-operation timeout is not a total acquisition deadline.

## Actual Isolated Execution

Export/copy `native-build` outside the checkout: 1,842 payloads,
43 scientific/1,798 workflow/one separately pinned setup overlay.
Payload bytes are 10,982,316; 1,821 Python files pass syntax checks.
The 745,425-byte source index SHA-256 is
`a65df12c5a0adfec4bbf5dc5b4e7b5dfa6913bb576fbc144d6ff40b161437682`.
Retain that digest independently before copied verification.

Four outer stages return zero: export, copied verification before, copied
acquisition, copied verification after. The copied stages use isolated
Python3.12, no inherited scientific environment and an empty PATH. No Git,
compiler, installed MAFFT/FastTree, source checkout or scientific dependency
is called by acquisition. Source modules and every acquired file are
rechecked; the source component stays unchanged.

All seven acquired files match the retained historical sizes/hashes. Every
file is mode 0644, including the FastTree ELF; no archive is extracted, binary
permission enabled, tool executed or native/scientific fixture rerun.
No existing installation, service, unrelated job or DGX is changed.
All artifacts remain local; no third-party source/binary is committed/uploaded.

The [machine receipt](publication_phylogeny_acquisition_20261002.json) is
25,745 bytes, SHA-256
`0e2d881ccede261b4bdd8d880742fed573fe450472c76fe1698b8796f34ffec0`.
It records the outer stage logs, anchored source, exact public downloads,
historical comparison and JUnit identity. The copied helper's 10,151-byte
completion receipt has SHA-256
`8f4b84bfe95c9fc5614db34c0c32ccbafdb5316fb6662d3fa5e6e219a76cf137`.

## Reader Command and Remaining Work

[Executable acquisition instructions](../PUBLICATION_SOURCE_COMPONENT.md#frozen-phylogeny-artifact-acquisition)
use fresh output outside an independently verified source component.
The next runtime step is a private MAFFT core build/staging controller that
does not require the historical installation or old fixture inputs. The
existing build/reconciliation admissions remain historical evidence and are
not automatically transferred to a different launcher or runtime.

The full MAFFT archive includes optional extension sources/notices even when
only core is later built. FastTree source/header and separate GPL license text
remain preserved; matching upstream binary bytes does not prove its compiler
recipe or redistribution compliance. Compiler/OS libraries, optional engines,
selected-file attribution/security review, full executor/handoff, remaining
QfO uncertainty, independent validation limitations, controlled resource
evidence and final manuscript/release/archive/DOI remain separate requirements.
Timing stays deferred; the overall publication goal is still active.
