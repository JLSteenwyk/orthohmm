"""Export the frozen canonical ordering comparison without replacing history."""

import argparse
import math
from pathlib import Path

from benchmark_tools.prepare_qfo_canonical_pairs import check_records
from benchmark_tools.run_qfo_order_replay import record, save
from benchmark_tools.run_simulation_methods import read_frozen

METRICS = ("GO", "EC", "VGNC", "SwissTrees", "TreeFam-A", "FAS")
CANONICAL_SHA = "194d38a5ef6448e9f2b81a1def944a51dd16885037f73e180b0330b4feb51cd6"
HISTORICAL_SHA = "26853f680d5858b09ed580d1bbf64dab0920d5eacd62bf84b4fb3ea988086145"


def compare(historical, canonical):
    if (historical["status"] != "historical_comparator_readmitted"
            or canonical["status"] != "canonical_qfo_assessment_admitted"
            or canonical["accuracy_admitted"] is not True):
        raise ValueError("Require independently admitted assessments")
    arms = [historical["assessment"], canonical["assessment"]]
    participants = ("ohmm_qfo_corrected_factorial_p1_c1_r1", "ohmm_qfo_canonical_20260927")
    for arm, participant in zip(arms, participants):
        if arm["participant"] != participant or set(arm["endpoints"]) != set(METRICS):
            raise ValueError("Wrong participant or endpoint inventory")
        scores = [arm["endpoints"][m]["score"] for m in METRICS]
        if any(type(s) not in (int, float) or not math.isfinite(s) or not 0 <= s <= 1 for s in scores):
            raise ValueError("Invalid endpoint score")
        if not math.isclose(sum(scores) / 6, arm["secondary_six_metric_mean"], rel_tol=0, abs_tol=1e-14):
            raise ValueError("Wrong secondary mean")
    rows = []
    for metric in METRICS:
        left, right = [a["endpoints"][metric] for a in arms]
        if left["axes"] != right["axes"] or left["score_semantics"] != right["score_semantics"]:
            raise ValueError("Changed endpoint definition")
        rows.append(dict(endpoint=metric, historical=left["score"], canonical=right["score"],
                         difference=right["score"] - left["score"],
                         semantics=left["score_semantics"], axes=left["axes"],
                         historical_native=left["native_participant"],
                         canonical_native=right["native_participant"],
                         sampling_confounded=metric == "FAS"))
    return rows


def export(repo, output):
    if output.exists():
        raise FileExistsError(output)
    canonical_path = repo / "benchmarks/work/qfo_canonical_assessment_20260927/admission.json"
    historical_path = repo / "benchmark_tools/results/qfo_canonical_historical_readmission_20260927.json"
    canonical = read_frozen(canonical_path, CANONICAL_SHA)
    historical = read_frozen(historical_path, HISTORICAL_SHA)
    checked = [record(canonical_path), record(historical_path), historical["original"], historical["fresh"],
               *canonical["checked_records"]]
    check_records(checked)
    rows = compare(historical, canonical)
    result = dict(status="canonical_qfo_ordering_comparison", rows=rows,
                  canonical_admission=record(canonical_path), historical_readmission=record(historical_path),
                  scheduler=canonical["scheduler"], native_tasks=len(canonical["native_tasks"]),
                  historical_mean=historical["assessment"]["secondary_six_metric_mean"],
                  canonical_mean=canonical["assessment"]["secondary_six_metric_mean"],
                  source=record(__file__), publication_ready=False,
                  limitations=canonical["limitations"])
    output.mkdir(parents=True)
    lines = ["# Canonical QfO Ordering Comparison", "",
             "| Endpoint | Historical retained | Canonical | Canonical minus historical |",
             "| --- | ---: | ---: | ---: |"]
    lines.extend(f"| {r['endpoint']} | {r['historical']:.9f} | {r['canonical']:.9f} | {r['difference']:+.9f} |" for r in rows)
    lines += ["", f"Secondary six-metric means: {result['historical_mean']:.9f} historical; "
              f"{result['canonical_mean']:.9f} canonical. These are not F1 scores.", "",
              "FAS uses unseeded native sampling; its difference is not an isolated ordering effect.",
              "GO/EC are mean Schlicker scores; VGNC/SwissTrees/TreeFam-A use harmonic TPR/PPV.",
              "Native axes, precision/recall and assessed-relation counts are retained in results.json.",
              "NR_ORTHOLOGS is not necessarily total submitted-pair coverage.",
              "No new paired confidence interval, independent generalization or default promotion.", ""]
    table = output / "scores.md"
    table.write_text("\n".join(lines))
    result["table"] = record(table)
    check_records(checked)
    save(output / "results.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    export(args.repo.resolve(), args.output.resolve())
