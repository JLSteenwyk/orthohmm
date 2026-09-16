"""Reproduce native normalization on fixed-length simulation search results.

Run with the installed OrthoFinder interpreter. This does not alter the tool,
input sequences, saved predictions, or scientific scores.
"""

import argparse
import csv
import gzip
import importlib.metadata
import json
from pathlib import Path
import sys
import warnings

import numpy as np
from scipy import sparse
from orthofinder.tools import waterfall

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.benchmark_production import file_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blast", type=Path, required=True)
    parser.add_argument("--query-count", type=int, required=True)
    parser.add_argument("--target-count", type=int, required=True)
    parser.add_argument("--length", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    if importlib.metadata.version("orthofinder") != "3.1.5" or min(args.query_count, args.target_count, args.length) < 1:
        raise ValueError("Expected OrthoFinder 3.1.5 and positive dimensions")
    matrix = sparse.lil_matrix((args.query_count, args.target_count))
    seen = set()
    with gzip.open(args.blast, "rt") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) != 12:
                raise ValueError("Expected native twelve-column DIAMOND output")
            i, j = int(row[0].split("_")[1]), int(row[1].split("_")[1])
            if (i, j) in seen:
                raise ValueError("Duplicate pair requires native aggregation review")
            seen.add((i, j))
            matrix[i, j] = float(row[11])
    lengths = {0: np.full(args.query_count, float(args.length)), 1: np.full(args.target_count, float(args.length))}
    li, lj, scores = waterfall.scnorm.GetLengthArraysForMatrix(matrix, lengths[0], lengths[1])
    top_l, top_s = waterfall.scnorm.GetTopPercentileOfScores(li * lj, scores, 95)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        parameters = waterfall.scnorm.CalculateFittingParameters(top_l, top_s)
        normalized = waterfall.WaterfallMethod.NormaliseScores(matrix, lengths, 0, 1).tocsr()
    source = Path(waterfall.__file__)
    result = {"schema_version": 1, "diagnostic_only": True, "accuracy_computed": False,
              "command": [sys.executable, *sys.argv], "python": sys.version,
              "packages": {n: importlib.metadata.version(n) for n in ("orthofinder", "numpy", "scipy")},
              "input": file_record(args.blast, args.blast.parent),
              "native_source": dict(file_record(source, source.parent), absolute_path=str(source)),
              "diagnostic_source": file_record(Path(__file__), Path(__file__).parent),
              "constant_length_assumption": args.length, "matrix_shape": list(matrix.shape),
              "raw_hits": matrix.nnz, "unique_length_products": np.unique(li * lj).tolist(),
              "fitting_parameters": parameters.tolist(), "normalized_stored": normalized.nnz,
              "nonfinite_normalized": int(np.count_nonzero(~np.isfinite(normalized.data))),
              "warnings": [str(w.message) for w in caught]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        handle.write(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
