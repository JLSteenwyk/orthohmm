# OrthoHMM Publication Handoff Candidate

This local candidate connects the retained manuscript review, frozen scientific
source, native-preparation workflow, current benchmark table and two standalone
arithmetic replays. It is not a submission-ready or redistribution-cleared
public release. Controlled timing and final runtime, rights and manuscript
reconciliation remain open. Scientific scores and settings are unchanged.

## Contents

| Location | Scope |
| --- | --- |
| `manuscript/` | Previously inspected nine-page main review, bibliography, direct local targets and review receipts |
| `source/scientific/` | Frozen scientific revision `7f3a9e40dd7e79f842cc2c11fb8b548f9a802806` |
| `source/workflow/` | Separately committed workflow/tests and four pinned acquisition/runtime support documents |
| `comparison/` | Current eight-tool table, 72 numeric cells and retained provenance |
| `arithmetic/` | Identifier-free YGOB counts, two simulation count reports, standalone scripts, guides and NumPy pin |

The manuscript is a dated snapshot, not a rendered copy of every subsequent
repository update. Open its relative HTML entrypoint from `REVIEW_INDEX.json`
for direct local links. The preserved PDF may retain historical path
annotations. Linked documents' own transitive assets are not all included.
The component guides explain their individual scope; their repository links
may refer to omitted evidence. Historical metadata paths are not rewritten.

## Build And Verify

From a full clone containing the frozen scientific/review revisions, select
an explicit committed workflow revision and a fresh output directory:

```bash
python -B benchmark_tools/bundle_publication_handoff.py build \
  --repo . --revision WORKFLOW_COMMIT --output /absolute/fresh/handoff
```

Keep the emitted `HANDOFF_INDEX.json` SHA-256 outside the directory before
transfer. After relocation, only standard-library Python is needed for integrity
and component verification, not Git, the original checkout or native tools:

```bash
python3 -I -S -B /relocated/handoff/bundle_publication_handoff.py verify \
  /relocated/handoff --manifest-sha256 RETAINED_HANDOFF_INDEX_SHA256
```

Every byte/mode and the exact inventory are checked before loading the included
standard-library component verifiers. An internal digest is not an external
authenticity anchor. Successful verification does not execute numerical replay,
re-admit native inference, validate all citations/anchors or establish rights.

## Numerical Replay

Use the recorded Python 3.12.3 and NumPy 2.2.6 in a private environment. NumPy
is the sole nonstandard import, not a bundled/hermetic compiled runtime. From
the handoff root, choose new environment/output paths outside this immutable
directory; its integrity verifier rejects extra generated files:

```bash
python3 -m venv /absolute/fresh/replay-env
/absolute/fresh/replay-env/bin/python -m pip install \
  -r arithmetic/requirements-ygob-arithmetic.txt
/absolute/fresh/replay-env/bin/python -I -B arithmetic/reproduce_ygob_validation.py reproduce \
  --snapshot arithmetic/ygob_sufficient_counts_20261002.json.gz \
  --snapshot-sha256 62885bb85eb357d97e55d053134d262e1ebb1e2b088a37df96986f98f8620d66 \
  --summary arithmetic/ygob_frozen_results_20260916.json --output /absolute/fresh/ygob-replay.json
/absolute/fresh/replay-env/bin/python -I -B arithmetic/reproduce_simulation_panels.py \
  --results arithmetic/simulation_fixed_native_results_20260916.json \
  --results arithmetic/simulation_variable_native_results_20260916.json \
  --output /absolute/fresh/simulation-replay.json
```

The parent directories must exist and the two output files must not. The
[YGOB guide](arithmetic/YGOB_ARITHMETIC_REPLAY.md) and
[simulation guide](arithmetic/SIMULATION_ARITHMETIC_REPLAY.md) retain statistic,
failure, bootstrap and development-exposure boundaries. These checks replay
retained-count arithmetic, not raw scoring, simulation truth, reference
independence or native inference. Fixed-length OrthoFinder comparisons stay
unavailable; failures are not imputed. The six-metric QfO mean is a project-
defined secondary summary, not official F1 or a universal tool ranking.

## Native Reproduction

Follow `source/workflow/benchmark_tools/PUBLICATION_SOURCE_COMPONENT.md` for
frozen OrthoBench acquisition verification, path rebinding and offline base
installation. Its controller requires separately supplied exact package archives,
pip bootstrap wheel and a trusted compatible Conda entrypoint. Historical
runtime reconstruction requires explicit acknowledgement; it is not a patched
general-purpose installation recommendation.

The integrated native workflow still needs separately provisioned scientific
wheel sets, aligner/tree builder, reader export and locks, and acquired raw
inputs. Its references are scoring inputs, not inference inputs. These are
not delivered by this source/count handoff. Do not execute historical DGX
scripts: the approved timing host is the local Threadripper, and production
timing requires its separate environment/resource/isolation gates.

No archive upload, public release, DOI, journal submission, raw-data rights
clearance, complete bootstrap/native/OS dependency closure or controlled
resource comparison is implied by this candidate.
