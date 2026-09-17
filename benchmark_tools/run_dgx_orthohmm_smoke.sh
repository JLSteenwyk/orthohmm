#!/usr/bin/env bash
# First frozen simulation dataset; functional portability, not timing evidence.
set -euo pipefail
root=$(realpath "${1:?Usage: run_dgx_orthohmm_smoke.sh PROJECT_ROOT}")
input="$root/smoke_missing20_20261101_input"
output="$root/orthohmm_pipeline_smoke_v4"
python="$root/envs/orthohmm/bin/python"
test "$(uname -m)" = aarch64
test ! -e "$output"
test -d "$input"
"$python" -c 'import importlib.metadata as m; assert m.version("DendroPy") == "5.0.8"'
cd "$root/core_arm_v2"
test "$(git rev-parse HEAD)" = 7f3a9e40dd7e79f842cc2c11fb8b548f9a802806
git diff --exit-code HEAD -- orthohmm
sha256sum --check <<'CHECKSUMS'
857a3832f9fb3a073257a499c700fcb16436db5d4917b20d67e4ad0cbc60c29f  orthohmm/search/csrc/hmm_viterbi.so
a0541ba3205217c60064965b1583b2d4a244096ff04a75c62891314549856973  orthohmm/search/csrc/kmer_prefilter.so
c7eaa03ebd3398f092c063d7b1869fbdf2fa0f1f478e320b4ca282823a2712c0  orthohmm/search/csrc/pair_align.so
CHECKSUMS
mkdir "$output"
mkdir "$output/high" "$output/satellite"
export PATH="$root/external-prefix-v1/mafft-7.525/bin:$root/external-prefix-v1/fasttree-2.2.0:$PATH"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
unset PYTHONPATH
sha256sum "$input"/*.fasta | tee "$output/input.before.sha256"
git ls-files -z orthohmm | xargs -0 sha256sum | tee "$output/source.before.sha256"
sha256sum orthohmm/search/csrc/*.so > "$output/native.before.sha256"
"$python" -m pip freeze > "$output/packages.txt"
settings=(-c 4 --threads_per_worker 4 -x BLOSUM62 -e 0.0001 --clustering leiden
          --cpm_resolution 0.1 --refinement_profile default --accuracy_profile high_sensitivity --stop infer)
set -x
high_status=0
"$python" -m orthohmm "$input" -o "$output/high" "${settings[@]}" \
    --metrics_json "$output/high.json" > "$output/high.log" 2>&1 || high_status=$?
satellite_status=0
"$python" -m orthohmm "$input" -o "$output/satellite" "${settings[@]}" \
    --metrics_json "$output/satellite.json" --phylogeny reconcile --species_tree_mode infer \
    --aligner mafft --tree_builder FastTree --phylogeny_candidates satellite_v2 \
    --phylogeny_root_rule species_overlap --phylogeny_pair_rule positive_paralogy \
    --species_tree_rooting min_variance > "$output/satellite.log" 2>&1 || satellite_status=$?
set +x
printf 'high_exit=%s\nsatellite_exit=%s\n' "$high_status" "$satellite_status" | tee "$output/exits.txt"
sha256sum --check "$output/input.before.sha256"
sha256sum --check "$output/source.before.sha256"
sha256sum --check "$output/native.before.sha256"
test "$high_status" = 0
test "$satellite_status" = 0
test -s "$output/high.json"
test -s "$output/satellite.json"
test -s "$output/high/orthohmm_orthogroups.txt"
test -s "$output/satellite/orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv"
