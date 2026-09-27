# Fixed Leiden-Distribution Downstream Replay

## Objective and Design

The private distribution control isolated an initial-clustering difference,
not its downstream effect. Run two fixed, sequential arms through the existing
saved-hit replay: clean leidenalg 0.12.0 and private 0.11.0 distribution. Both
use the fresh installed OrthoBench numeric search checkpoint and its same 12
proteomes (251,378 genes). No search is repeated; no reference-label path is
passed to inference. No parameters, seeds, profiles or endpoints are selected
from the results.

Hold all scientific source and other installed packages fixed. Use CPM 0.1,
seed 4, BLOSUM62, one profile iteration, minimum one species per profile, 32
CPUs and one BLAS/OpenMP thread. Preserve initial, multipass, profile-base and
profile-final clustering graph/partition snapshots. Retain multipass and
profile refinement outputs, then apply the frozen satellite_v2 candidate
expansion. Stop before phylogeny and reference scoring.

Private venvs use path-only `.pth` files pointing to the unchanged clean
site-packages and, for the 0.11 arm, the validated private distribution first.
This is a diagnostic environment, not a portable installation or a proposed
dependency downgrade. Parent scientific imports are loaded before exposing
audit scripts. A plain-child import probe must resolve the same scientific
module and native Leiden extension; this matters because production clustering
spawns a new Python worker. Paths, sources and input hashes are frozen in the
execution plan. Do not edit the driver or replay helper while this job runs.

Each arm has one attempt and a six-hour timeout; a failed arm stops the sequence.
Submit with 32 CPUs, 128 GiB RAM and 13 hours total. No resume or automatic
retry. Preserve logs and failed outputs. This is shared-host diagnostic work,
not the dedicated matched-resource scaling panel. DGX remains unused.

## Preparation and Validation

The [preparation receipt](ob_dependency_replay_preparation_20260926.json)
records the pinned v2 plan and two passing small complete saved-hit replays,
including four clustering calls and candidate expansion in each distribution.
The child import probe initially rejected a smoke call launched from repository
cwd, before clustering: its plain child resolved worktree source. The failure
is preserved; corrected fixtures use the same private cwd as the production
controller. The earlier preparation directory/plan was never submitted and
remains retained. No full-dataset retry occurred.

V2 plan:
`benchmarks/work/ob_dependency_replay_v2_20260926/plan.json`

SHA-256:
`950a472f1ca767851c69446412edc0b26ea5fc5646160c1d385f6318327b87ba`

```sh
python -m benchmark_tools.run_ob_dependency_replay --execute \
  benchmarks/work/ob_dependency_replay_v2_20260926/plan.json \
  --plan-sha256 950a472f1ca767851c69446412edc0b26ea5fc5646160c1d385f6318327b87ba
```

After successful native completion, independently verify both execution records,
input/source hashes, package identities and exactly four saved clustering calls.
Check every partition for unique complete coverage; compare gene memberships
without depending on group labels. Compare candidate partitions with each other
and retained historical/fresh candidates. Any later score readback must retain
both arms, be labeled post-hoc/development-exposed, and must not replace current
benchmark tables or support a new superiority claim. Causal attribution to a
final phylogenetic F1 change requires additional evidence beyond these pre-tree
results. The broader publication goal remains incomplete.
