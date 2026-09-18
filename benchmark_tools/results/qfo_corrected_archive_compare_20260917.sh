#!/bin/bash
#SBATCH --job-name=qfo_corrected_compare
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=benchmarks/work/qfo_corrected_source_20260917/compare_%j.log
set -euo pipefail

root=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
executor="$root/benchmarks/work/publication_qfo_corrected_archive_audit_v1"
revision=a1fbffcc8cb5b4fdabca18843d79d94b0baadaef
archive="$root/benchmarks/work/qfo_corrected_source_20260917/QfO_release_2020_04_with_updated_UP000008143.tar.gz"
output="$root/benchmark_tools/results/qfo_corrected_archive_comparison_20260917.json"

test "$(git -C "$executor" rev-parse HEAD)" = "$revision"
git -C "$executor" diff --exit-code HEAD -- benchmark_tools
test "$(stat -c %s "$archive")" -eq 2648666198
test ! -e "$output"
accounting=$(sacct -j 21687 --parsable2 --noheader --format=JobID,State,ExitCode,NodeList)
printf '%s\n' "$accounting"
test "$(printf '%s\n' "$accounting" | awk -F '|' '$1 == "21687" {print $2 "|" $3 "|" $4}')" = 'COMPLETED|0:0|bizon'
printf 'Comparison job %s starts %s; executor %s\n' "$SLURM_JOB_ID" "$(date -u +%FT%TZ)" "$revision"
export PYTHONHASHSEED=0 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$root"
/home/bizon/anaconda3/bin/python -c 'import sys, Bio; print(sys.version); print("Biopython", Bio.__version__)'
/home/bizon/anaconda3/bin/python "$executor/benchmark_tools/compare_qfo_corrected_archive.py" \
  --archive "$archive" \
  --prepared "$root/benchmark_tools/results/qfo_factorial_prepared_20260917.json" \
  --mapping "$root/qfo_benchmark/benchmark-webservice/reference_data/2020/mapping.json.gz" \
  --aliases "$root/benchmark_tools/results/swiss_sequence_alias_audit_20260917.json" \
  --output "$output"
sha256sum "$output"
printf 'Comparison completed %s; independent interpretation still required.\n' "$(date -u +%FT%TZ)"
