"""Compare GO/EC scored-pair sets for all eight frozen corrected methods."""

import argparse
from itertools import combinations
import json
from pathlib import Path

from benchmark_tools.compare_qfo_scored_pairs import compare, read_scores
from benchmark_tools.map_corrected_vgnc_blocks import MANIFEST, MANIFEST_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def run(repo, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    identity = record(repo / MANIFEST)
    if identity["sha256"] != MANIFEST_SHA:
        raise ValueError("Changed corrected comparison manifest")
    manifest = json.loads((repo / MANIFEST).read_text())
    methods = manifest["methods"]
    if (len(methods) != 8 or len({m["key"] for m in methods}) != 8
            or any(m["status"] != "admitted" for m in methods)):
        raise ValueError("Require eight distinct admitted methods")
    checked = [identity, record(__file__), record(Path(__file__).with_name("compare_qfo_scored_pairs.py"))]
    inputs = {}
    for method in methods:
        check(method["admission"])
        checked.append(method["admission"])
        admission = json.loads(Path(method["admission"]["path"]).read_text())
        for metric in ("GO", "EC"):
            candidates = [r for r in admission["metric_files"]
                          if r["path"].endswith(f"/results/{metric}/{metric}.json")]
            if len(candidates) != 1:
                raise ValueError("Ambiguous admitted endpoint")
            endpoint = candidates[0]
            check(endpoint)
            checked.append(endpoint)
            paths = list(Path(endpoint["path"]).parent.glob("*raw.txt.gz"))
            if len(paths) != 1:
                raise ValueError("Ambiguous raw scored-pair file")
            raw = record(paths[0])
            checked.append(raw)
            inputs[method["key"], metric] = raw
    rows, summaries = [], []
    for metric in ("GO", "EC"):
        scores = {}
        for method in methods:
            key = method["key"]
            values = read_scores(Path(inputs[key, metric]["path"]), metric)
            mean = sum(values.values()) / (len(values) * 1e6)
            if (len(values) != method["details"][metric]["assessed_relations"]
                    or abs(mean - method["scores"][metric]) > 5.051e-7):
                raise ValueError("Raw count/mean differs from admitted endpoint")
            scores[key] = values
            summaries.append(dict(metric=metric, method=key, raw=inputs[key, metric],
                assessed_pairs=len(values), rounded_mean=mean, admitted_mean=method["scores"][metric]))
        for left, right in combinations(methods, 2):
            rows.append(dict(metric=metric, left=left["key"], right=right["key"],
                result=compare(scores[left["key"]], scores[right["key"]])))
        del scores
    for item in checked:
        check(item)
    result = dict(status="corrected_all_method_scored_pair_panel", checked_records=checked,
        methods=[dict(key=m["key"], label=m["label"], prediction_semantics=m["prediction_semantics"])
                 for m in methods], endpoints=summaries, comparisons=rows,
        uncertainty_admitted=False, publication_ready=False,
        limitations=["Descriptive comparison of eligible scored pairs, not all predictions.",
            "Identical serialized scores do not prove identical unavailable full-precision scores.",
            "Intersection-only means are not replacements for native endpoint means.",
            "No family independence, paired confidence interval or causal attribution established."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.repo.resolve(), args.output)
