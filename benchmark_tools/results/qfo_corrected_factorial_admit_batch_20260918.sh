#!/bin/bash
#SBATCH --job-name=qfo_corrected_factorial_admit
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --output=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/qfo_corrected_factorial_admit_%j.log
set -euo pipefail
ROOT=/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm
EXECUTOR=${1:?Frozen admission executor required}
COMMIT=${2:?Exact revision required}
INDEX=${3:?Reconciliation index required}
ARRAY=${4:?Reconciliation array job required}
ADMISSION_JOB=${5:?Candidate admission job required}
[[ $# == 5 && "$COMMIT" =~ ^[0-9a-f]{40}$ && "$INDEX" =~ ^[0-3]$ ]]
[[ "$ARRAY" =~ ^[0-9]+$ && "$ADMISSION_JOB" =~ ^[0-9]+$ ]]
[[ $(git -C "$EXECUTOR" rev-parse HEAD) == "$COMMIT" ]]
git -C "$EXECUTOR" diff --quiet HEAD -- benchmark_tools orthohmm
JOB=$(/home/bizon/anaconda3/bin/python - "$ARRAY" "$INDEX" <<'PY'
import csv
import io
import subprocess
import sys

task = sys.argv[1] + "_" + sys.argv[2]
accounting = subprocess.check_output([
    "sacct", "-j", task, "--parsable2",
    "--format=JobID%64,JobIDRaw%64,State,ExitCode"], text=True)
rows = [row for row in csv.DictReader(io.StringIO(accounting), delimiter="|")
        if row["JobID"] == task]
if len(rows) != 1:
    raise ValueError("Require one exact array-task accounting record")
row = rows[0]
if (row["State"], row["ExitCode"]) != ("COMPLETED", "0:0"):
    raise ValueError("Reconciliation task did not complete successfully")
if not row["JobIDRaw"].isascii() or not row["JobIDRaw"].isdigit():
    raise ValueError("Expected numeric raw task job ID")
print(row["JobIDRaw"])
PY
)
ADMISSION="$ROOT/benchmarks/work/qfo_corrected_candidate_admission_20260918.json"
read -r SHA _ < <(sha256sum "$ADMISSION")
OUTPUT="$ROOT/benchmarks/work/qfo_corrected_factorial_native_admission_${INDEX}_20260918.json"
export PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
exec /home/bizon/anaconda3/bin/python "$EXECUTOR/benchmark_tools/admit_qfo_corrected_factorial_cell.py" \
    --root "$ROOT" --index "$INDEX" --job "$JOB" --admission "$ADMISSION" \
    --admission-sha256 "$SHA" --admission-job "$ADMISSION_JOB" --output "$OUTPUT"
