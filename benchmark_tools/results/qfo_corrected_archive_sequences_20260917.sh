#!/bin/bash
#SBATCH --job-name=qfo_corrected_sequences
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=benchmarks/work/qfo_corrected_source_20260917/sequences_%j.log
set -euo pipefail
root=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
executor="$root/benchmarks/work/publication_qfo_archive_sequences_v1"
test "$(git -C "$executor" rev-parse HEAD)" = 25394c3e2b5b7edf92b710d5edb6859b1d09866c
git -C "$executor" diff --exit-code HEAD -- benchmark_tools
archive="$root/benchmarks/work/qfo_corrected_source_20260917/QfO_release_2020_04_with_updated_UP000008143.tar.gz"
test "$(stat -c %s "$archive")" -eq 2648666198
accounting=$(sacct -j 21687,21688 --parsable2 --noheader --format=JobID,State,ExitCode,NodeList)
printf '%s\n' "$accounting"
for job in 21687 21688; do
  test "$(printf '%s\n' "$accounting" | awk -F '|' -v job="$job" '$1 == job {print $2 "|" $3 "|" $4}')" = 'COMPLETED|0:0|bizon'
done
export PYTHONHASHSEED=0 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$root"
printf 'Sequence compatibility audit %s starts %s\n' "$SLURM_JOB_ID" "$(date -u +%FT%TZ)"
/home/bizon/anaconda3/bin/python -c 'import sys, Bio; print(sys.version); print("Biopython", Bio.__version__)'
/home/bizon/anaconda3/bin/python "$executor/benchmark_tools/audit_qfo_archive_sequences.py" \
  --archive "$archive" \
  --prepared "$root/benchmark_tools/results/qfo_factorial_prepared_20260917.json" \
  --mapping "$root/qfo_benchmark/benchmark-webservice/reference_data/2020/mapping.json.gz" \
  --database "$root/qfo_benchmark/benchmark-webservice/reference_data/2020/ServerIndexed.db" \
  --output "$root/benchmark_tools/results/qfo_corrected_archive_sequences_20260917.json"
sha256sum "$root/benchmark_tools/results/qfo_corrected_archive_sequences_20260917.json"
printf 'Audit finished %s; report requires independent interpretation.\n' "$(date -u +%FT%TZ)"
