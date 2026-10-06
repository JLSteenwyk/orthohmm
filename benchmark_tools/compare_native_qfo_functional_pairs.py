"""Describe scored-pair composition for the two admitted native P0/C0 QfO cells."""

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import statistics
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import export_native_qfo_scientific_scores as reporter
from benchmark_tools.audit_qfo_fas_samples import read_sample
from benchmark_tools.compare_qfo_scored_pairs import compare, read_scores
from benchmark_tools.export_native_factorial_progress import load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

CELLS = ("p0_c0_r1", "p0_c0_r0")


def raw_bindings(row, evidence):
    admission, admission_ref = load(row["admission"]["path"], row["admission"]["sha256"], evidence)
    require(admission_ref == row["admission"] and admission["accuracy_admitted"] is True
            and admission["cell"] == row["cell"] and admission["native_index"] == row["index"]
            and admission["native_job_id"] == row["native_job_id"], "Changed native admission identity")
    execution_ref = admission["execution_report"]
    require(execution_ref in admission["checked_records"], "Execution absent from admission records")
    execution, observed = load(execution_ref["path"], execution_ref["sha256"], evidence)
    require(observed == execution_ref and execution["exit_code"] == 0
            and execution["status"] == "process_succeeded_pending_independent_admission",
            "Require successful independently admitted assessment")
    result = {}
    for metric in ("GO", "EC", "FAS"):
        paths = [ref for ref in execution["outputs"] if f"/results/{metric}/" in ref["path"]
                 and ref["path"].endswith("raw.txt.gz")]
        require(len(paths) == 1 and paths[0] in admission["checked_records"],
                "Ambiguous or unadmitted raw " + metric)
        check(paths[0])
        evidence.append(paths[0])
        result[metric] = paths[0]
    require(result["FAS"] == admission["fas_sample"]["raw"], "FAS sample identity differs")
    return admission, result


def protein_structure(pairs):
    degrees = Counter(protein for pair in pairs for protein in pair)
    return dict(proteins=len(degrees), proteins_in_multiple_pairs=sum(v > 1 for v in degrees.values()),
                maximum_protein_degree=max(degrees.values()))


def fas_overlap(left, right):
    shared = left.keys() & right.keys()
    differences = [left[p] - right[p] for p in shared]
    nl, nr, ns = len(left), len(right), len(shared)
    require(nl >= 2 and nr >= 2, "Insufficient FAS samples")
    left_sum, right_sum = math.fsum(left.values()), math.fsum(right.values())
    shared_left, shared_right = math.fsum(left[p] for p in shared), math.fsum(right[p] for p in shared)
    return dict(left_sample_pairs=nl, right_sample_pairs=nr, shared_sample_pairs=ns,
        left_only_sample_pairs=nl-ns, right_only_sample_pairs=nr-ns,
        shared_fraction_of_left=ns/nl, shared_fraction_of_right=ns/nr,
        shared_pairs_with_different_serialized_scores=sum(d != 0 for d in differences),
        maximum_shared_absolute_difference=max(map(abs, differences), default=None),
        shared_conditional_mean_difference=math.fsum(differences)/ns if ns else None,
        original_sample_mean_difference=left_sum/nl-right_sum/nr,
        original_sample_mean_difference_components=dict(
            shared_with_original_denominators=shared_left/nl-shared_right/nr,
            left_only=(left_sum-shared_left)/nl, negative_right_only=-(right_sum-shared_right)/nr))


def run(snapshot_path, snapshot_sha, output):
    require(not output.exists() and not output.is_symlink(), "Output already exists")
    evidence = []
    snapshot, snapshot_ref = load(snapshot_path, snapshot_sha, evidence)
    require(snapshot["schema"] == "native_qfo_scientific_reporting_snapshot_v1"
            and snapshot["source"] == record(reporter.__file__) and snapshot["publication_ready"] is False,
            "Require source-bound incomplete native scientific snapshot")
    normal, recovered = [], []
    for row in snapshot["rows"]:
        destination = (normal if row["status"] == "supplied_native_admission" else recovered
                       if row["status"] == "supplied_recovered_scientific_admission" else None)
        if destination is not None:
            destination.append((row["admission"]["path"], row["admission"]["sha256"]))
    replay = reporter.collect(snapshot["plan"]["path"], snapshot["plan"]["sha256"], normal, recovered)
    require(all(snapshot[k] == value for k, value in replay.items()), "Scientific snapshot replay differs")
    rows = {row["cell"]: row for row in snapshot["rows"]}
    require(len(rows) == len(snapshot["rows"]) and all(rows[c]["accuracy_admitted"] is True for c in CELLS),
            "Require two distinct admitted native P0/C0 cells")
    repaired = rows[CELLS[0]]
    require(repaired["status"] == "supplied_recovered_scientific_admission"
            and repaired["resources"] is None and repaired["timing_admitted"] is False
            and repaired["timing_eligible"] is False, "Recovered accuracy relabels failed timing")
    admissions, raws = {}, {}
    for cell in CELLS:
        admissions[cell], raws[cell] = raw_bindings(rows[cell], evidence)
    summaries, comparisons = [], []
    for metric in ("GO", "EC", "FAS"):
        scores = {}
        for cell in CELLS:
            raw = raws[cell][metric]
            values = (read_sample(Path(raw["path"])) if metric == "FAS" else read_scores(Path(raw["path"]), metric))
            scores[cell] = values
            mean = statistics.mean(values.values()) if metric == "FAS" else sum(values.values())/(len(values)*1e6)
            tolerance = 1e-12 if metric == "FAS" else 5.051e-7
            expected_count = (admissions[cell]["fas_sample"]["sample_pairs"] if metric == "FAS"
                              else rows[cell]["endpoint_details"][metric]["assessed_relations"])
            require(len(values) == expected_count and abs(mean-rows[cell]["scores"][metric]) <= tolerance,
                    "Raw count or mean differs from admitted " + metric)
            summaries.append(dict(cell=cell, metric=metric, raw=raw, scored_pairs=len(values),
                native_mean=rows[cell]["scores"][metric], raw_mean=mean, mean_tolerance=tolerance,
                **protein_structure(values)))
        result = (fas_overlap(scores[CELLS[0]], scores[CELLS[1]]) if metric == "FAS"
                  else compare(scores[CELLS[0]], scores[CELLS[1]]))
        comparisons.append(dict(metric=metric, left=CELLS[0], right=CELLS[1], result=result))
        del scores
    evidence.extend([record(reporter.__file__), record(__file__), *[record(Path(__file__).with_name(name))
                     for name in ("compare_qfo_scored_pairs.py", "audit_qfo_fas_samples.py")]])
    for ref in evidence:
        check(ref)
    result = dict(schema="native_qfo_functional_pair_composition_v1", snapshot=snapshot_ref,
        source=record(__file__), evidence=evidence, cells=list(CELLS), endpoints=summaries,
        comparisons=comparisons, new_scoring_or_admission=False, new_bootstrap_draws=0,
        uncertainty_admitted=False, scientific_timings_admitted=False, publication_ready=False,
        limitations=["Native P0/C0 two-cell diagnostic; initial HMM search on, not a selected-default comparison.",
            "GO/EC use six-decimal serialized eligible scores, not all predicted relations or annotation recomputation.",
            "FAS compares retained unseeded realized samples, not full eligible prediction populations.",
            "Shared-pair conditioning changes the endpoint; original means and denominators remain unchanged.",
            "Overlap or repeated-protein counts do not define independent units or validate paired intervals.",
            "Decomposition is arithmetic, not causal attribution, biological accuracy or general superiority.",
            "Recovered native timing remains failed/null/ineligible; no native rerun or uncertainty admission."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("snapshot", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--snapshot-sha256", required=True)
    args = parser.parse_args()
    report = run(args.snapshot, args.snapshot_sha256, args.output)
    print(json.dumps(dict(endpoints=len(report["endpoints"]), comparisons=len(report["comparisons"]),
                         uncertainty_admitted=report["uncertainty_admitted"])))
