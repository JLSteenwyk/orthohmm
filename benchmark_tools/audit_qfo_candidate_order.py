"""Read-only ordering compatibility audit of the admitted corrected QfO hits."""

import argparse
import json
from pathlib import Path

import numpy as np

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

ADMISSION_SHA = "7ae64b1ccd398a7c61c6011f89582c37be91fb9b3954edf85f0cf7d956890f2d"


def order_statistics(queries, targets, scores, genes, chunk_size=1000000):
    if chunk_size < 1 or genes < 1:
        raise ValueError("Positive chunk size and gene universe required")
    if any(x.ndim != 1 for x in (queries, targets, scores)) or not len(queries) == len(targets) == len(scores):
        raise ValueError("Invalid hit shapes")
    if queries.dtype.kind not in "iu" or targets.dtype.kind not in "iu" or scores.dtype.kind not in "fiu":
        raise ValueError("Invalid hit dtypes")
    descents = duplicates = self_hits = 0
    previous = None
    for start in range(0, len(queries), chunk_size):
        q, t, s = (x[start:start+chunk_size] for x in (queries, targets, scores))
        if np.any(q < 0) or np.any(t < 0) or np.any(q >= genes) or np.any(t >= genes):
            raise ValueError("Hit index outside universe")
        if not np.isfinite(s).all() or np.any(s <= 0):
            raise ValueError("Scores must be finite and positive")
        descents += int(np.count_nonzero((q[1:] < q[:-1]) | ((q[1:] == q[:-1]) & (t[1:] < t[:-1]))))
        duplicates += int(np.count_nonzero((q[1:] == q[:-1]) & (t[1:] == t[:-1])))
        first, last = (int(q[0]), int(t[0])), (int(q[-1]), int(t[-1]))
        if previous is not None:
            descents += first < previous
            duplicates += first == previous
        previous = last
        self_hits += int(np.count_nonzero(q == t))
    return dict(hits=len(queries), adjacent_descents=int(descents), adjacent_duplicates=int(duplicates),
        self_hits=self_hits, canonical_input_unchanged=descents == 0 and duplicates == 0,
        global_uniqueness_proven=descents == 0 and duplicates == 0)


def audit(repo):
    admission_path = repo / "benchmark_tools/results/qfo_corrected_candidate_admission_21759.json"
    admission_record = record(admission_path)
    if admission_record["sha256"] != ADMISSION_SHA:
        raise ValueError("Changed historical admission")
    admission = json.loads(admission_path.read_text())
    if admission["status"] != "corrected_qfo_candidates_admitted":
        raise ValueError("Missing admitted corrected candidates")
    prepared_record = admission["prepared_manifest"]
    check(prepared_record)
    prepared = json.loads(Path(prepared_record["path"]).read_text())
    checkpoint = prepared["numeric_checkpoint"]["manifest"]
    check(checkpoint)
    path = Path(checkpoint["path"])
    manifest = json.loads(path.read_text())
    if manifest["complete"] is not True or manifest["genes"] != 984137:
        raise ValueError("Wrong corrected checkpoint")
    records = [admission_record, prepared_record, checkpoint]
    for name, identity in manifest["files"].items():
        if Path(name).name != name:
            raise ValueError("Invalid checkpoint filename")
        item = dict(path=str(path.parent / name), **identity)
        check(item)
        records.append(item)
    names = (path.parent / "gene_names.txt").read_text().splitlines()
    if len(names) != manifest["genes"] or any(a >= b for a, b in zip(names, names[1:])):
        raise ValueError("Require unique lexically sorted gene names")
    arrays = [np.load(path.parent / name, mmap_mode="r", allow_pickle=False)
              for name in ("hit_queries.npy", "hit_targets.npy", "hit_scores.npy")]
    result = order_statistics(*arrays, len(names))
    if result["hits"] != manifest["hits"]:
        raise ValueError("Hit count differs")
    for item in records:
        check(item)
    return dict(status="admitted_qfo_candidate_input_order_audited", source=record(__file__),
        checked_records=records, genes=len(names), result=result, publication_ready=False,
        accuracy_evaluated=False, native_inference_executed=False,
        limitations=["Only ordering/array prerequisites at the retained candidate input are checked.",
            "No equivalence of installed dependencies, upstream inference, candidate outputs or phylogeny follows.",
            "If descents exist, adjacent duplicate counts do not prove global uniqueness."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.repo.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(result["result"], indent=2))
