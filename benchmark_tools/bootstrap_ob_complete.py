"""Exploratory paired uncertainty for all eight retained OrthoBench methods."""

import argparse
import json
from pathlib import Path

from benchmark_tools.bootstrap_orthobench import paired_bootstrap
from benchmark_tools.export_ob_complete_strata import PINS
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

BASELINE = "orthofinder_3_1_5_full"


def load_scores(base):
    names = ("orthobench_paired_uncertainty_20260916.json", "retained_ob_comparator_readback_20260926.json")
    refs, data = [], []
    for name in names:
        ref = record(base / name)
        if ref["sha256"] != PINS[name]:
            raise ValueError("Changed retained result")
        refs.append(ref)
        data.append(json.loads((base / name).read_text()))
    main, others = data
    if not others["all_scores_agree"]:
        raise ValueError("Comparator scores not verified")
    scores = dict(main["scores"])
    for row in others["rows"]:
        if row["key"] in scores or row["agrees_with_retained_score"] is not True:
            raise ValueError("Duplicate or unverified method")
        scores[row["key"]] = row["score"]
    if len(scores) != 8:
        raise ValueError("Require complete panel")
    refs.extend(others["checked_records"])
    refs.extend([*main["inputs"]["predictions"].values(), *main["inputs"]["references"], *main["inputs"]["uncertain"]])
    for ref in refs:
        check(ref)
    return scores, refs


def run(base, output, protocol_sha):
    if output.exists():
        raise FileExistsError(output)
    protocol = record(base / "OB_COMPLETE_UNCERTAINTY_PROTOCOL_20260928.md")
    if protocol["sha256"] != protocol_sha:
        raise ValueError("Protocol changed")
    scores, refs = load_scores(base)
    refs.extend([protocol, record(Path(__file__).resolve()),
                 record(Path(__file__).with_name("bootstrap_orthobench.py").resolve())])
    result = paired_bootstrap({m: s["refog_records"] for m, s in scores.items()}, BASELINE,
                              replicates=100000, seed=20260928, multiplicity_endpoints=21)
    if len(result["families"]) != 70:
        raise ValueError("Wrong reference panel")
    for method, point in result["point_estimates_percent"].items():
        if any(abs(value - scores[method][metric]) > 1e-8 for metric, value in point.items()):
            raise ValueError("Point estimates differ from retained scores")
    for ref in refs:
        check(ref)
    result.update(status="complete_ob_exploratory_paired_uncertainty", checked_records=refs,
                  publication_ready=False, independent_confirmation=False)
    with output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    args = parser.parse_args()
    result = run(args.base.resolve(), args.output, args.protocol_sha256)
    print(json.dumps(dict(methods=len(result["point_estimates_percent"]), contrasts=len(result["comparisons"]))))
