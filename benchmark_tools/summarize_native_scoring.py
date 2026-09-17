"""Describe retained portability discrepancies without relaxing the admission gate."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.probe_native_scoring import BANDS, compare, digest, fixtures


def summarize(left, right):
    queries, targets, pairs = fixtures()
    fixture = {"queries": [q.tolist() for q in queries], "targets": [t.tolist() for t in targets],
               "pairs": pairs.tolist(), "bands": list(BANDS)}
    for data in (left, right):
        if data["fixture_sha256"] != digest(fixture) or data["pairs"] != len(pairs):
            raise ValueError("Unexpected diagnostic fixtures")
        if [r["band"] for r in data["bands"]] != list(BANDS):
            raise ValueError("Incomplete diagnostic bands")
        for row in data["bands"]:
            for key in ("scores", "scalar_scores", "jit_scores", "evalues", "normalized"):
                if len(row[key]) != len(pairs):
                    raise ValueError("Incomplete diagnostic scores")
    for key in ("probe_sha256", "commit", "packages"):
        if left[key] != right[key]:
            raise ValueError("Diagnostic provenance differs")
    result = {"status": "diagnostic_only", "admitted": False, "pairs_per_band": len(pairs),
              "backends": [left["backend"], right["backend"]], "bands": []}
    try:
        result["gate_result"] = compare(left, right)
    except ValueError as error:
        result["gate_error"] = str(error)
    for a, b in zip(left["bands"], right["bands"]):
        differences = []
        for i, (x, y) in enumerate(zip(a["scores"], b["scores"])):
            if x != y:
                q, t = pairs[i]
                differences.append({"pair_index": i, "query_index": int(q), "target_index": int(t),
                                    "query_length": len(queries[q]), "target_length": len(targets[t]),
                                    "left_score": x, "right_score": y})
        result["bands"].append({"band": a["band"], "native_score_differences": differences,
                                "scalar_scores_exact": a["scalar_scores"] == b["scalar_scores"],
                                "jit_scores_exact": a["jit_scores"] == b["jit_scores"],
                                "evalues_exact": a["evalues"] == b["evalues"],
                                "normalized_exact": a["normalized"] == b["normalized"],
                                "decisions_exact": a["decisions"] == b["decisions"]})
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs=2, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    blobs = [p.read_bytes() for p in args.inputs]
    result = summarize(*(json.loads(blob) for blob in blobs))
    result["inputs"] = [{"path": str(p.resolve()), "sha256": hashlib.sha256(blob).hexdigest()}
                        for p, blob in zip(args.inputs, blobs)]
    with args.output.open("x") as handle:
        handle.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
