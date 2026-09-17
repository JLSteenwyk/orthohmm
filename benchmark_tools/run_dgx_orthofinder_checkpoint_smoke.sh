#!/usr/bin/env bash
# Historical filename: -og in 3.1.5 also runs phylogeny. Not sequence-only/timing.
set -euo pipefail
root=$(realpath "${1:?Usage: run_dgx_orthofinder_checkpoint_smoke.sh PROJECT_ROOT RECIPE_DIRECTORY}")
# Slurm executes a spooled script; $0 does not locate adjacent helpers.
recipe=$(realpath "${2:?Require original recipe directory}")
input="$root/smoke_missing20_20261101_input"
output="$root/orthofinder_checkpoint_smoke_v2"
prefix="$root/envs/orthofinder"
test "$(uname -m)" = aarch64
test ! -e "$output"
test -d "$input"
export PATH="$prefix/bin:$root/diamond-2.0.13-build-v2:$root/mcl-prefix-v3/bin:$root/external-prefix-v1/fasttree-2.1.11:$root/external-prefix-v1/mafft-7.525/bin:$root/famsa-source-v1:/usr/bin:/bin"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
unset PYTHONPATH CONDA_PREFIX
"$prefix/bin/python" -c 'import importlib.metadata as m; assert m.version("orthofinder") == "3.1.5"'
mkdir "$output"
mkdir "$output/input"
cp "$input"/*.fasta "$output/input/"
cd "$root/core_arm_v2"
sha256sum "$input"/*.fasta > "$output/input.before.sha256"
sha256sum "$output/input"/*.fasta > "$output/copies.before.sha256"
"$prefix/bin/python" "$recipe/inspect_orthofinder_runtime.py" --output "$output/runtime.before.json"
"$prefix/bin/python" - "$output/runtime.before.json" <<'PY'
import json
import sys
report = json.load(open(sys.argv[1]))
expected = {
    "diamond": "b79dfffa2878abcfa613f629e1b5da8e2994d2c76b6747f74b72ed6237758373",
    "FastTree": "8f9db34012ade374407809c275b909092b2ca1309726d5beaa52a0ad4294e3db",
    "famsa": "0a3988e9f7bae4dce39c730c0ceb3212d3ec9462ee585fe660f3861cbd38a7bb",
}
for name, digest in expected.items():
    if report["child_tools"][name]["file"]["sha256"] != digest:
        raise ValueError(f"Unexpected child executable: {name}")
mcl = report["child_tools"]["mcl"]
if mcl["exit_code"] != 0 or "14-137" not in mcl["stdout"] + mcl["stderr"]:
    raise ValueError("Unexpected MCL version")
PY
status=0
strace -f -s 16384 -e trace=execve -o "$output/execve.log" \
    "$prefix/bin/orthofinder" -f "$output/input" -t 4 -a 4 -S diamond \
    -og -o "$output/results" > "$output/native.log" 2>&1 || status=$?
printf 'exit=%s\n' "$status" | tee "$output/exit.txt"
sha256sum --check "$output/input.before.sha256"
sha256sum --check "$output/copies.before.sha256"
"$prefix/bin/python" "$recipe/inspect_orthofinder_runtime.py" --output "$output/runtime.after.json"
cmp "$output/runtime.before.json" "$output/runtime.after.json"
test "$status" = 0
test -s "$output/execve.log"
test "$(find "$output/results" -name 'clusters_OrthoFinder_I*.txt_id_pairs.txt' | wc -l)" = 1
test "$(find "$output/results" -name 'OrthoFinder_graph.txt' | wc -l)" = 1
