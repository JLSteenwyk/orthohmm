#!/bin/sh
# This script was created by sbatch --wrap.

exec env -u PYTHONHOME -u PYTHONPATH -u LD_PRELOAD -u LD_LIBRARY_PATH -u LD_AUDIT PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/native_factorial_review_py310_20261004/bin/python -B -m benchmark_tools.run_native_factorial_qfo_assessment --root /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm --pairs /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/native_factorial_qfo_pairs_22435/results.json --pairs-sha256 be4243d8470c9dd4e7f23a8d734824aab8b67c63a1acbd54c958868b9ec5009d --conversion-job 22439
