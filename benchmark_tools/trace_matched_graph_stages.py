"""Post-hoc complete-panel localization of matched-search graph pair errors."""

import argparse
import json
from pathlib import Path
from statistics import mean

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from benchmark_tools.audit_matched_graph import partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.score_matched_graph import CONDITIONS, SEEDS
from benchmark_tools.simulation_conditions import group_pairs, score_pairs


READBACK_SHA = "bbc9e6b40c75682e0336e3a466b8c2e740eed54f4b6a762203984bfa8c673a9f"
SCORES_SHA = "d31dbc0f4e3f02dfb534c8b59b559db00260227f23f2640b4fed13d5aee8ec89"
STAGES = ("direct_hits", "rbnh_edges", "rbnh_components", "initial",
          "multipass_edges", "multipass_components", "multipass", "final")


def edge_pairs(names, owners, sources, targets):
    if len(sources) != len(targets):
        raise ValueError("Unequal edge arrays")
    result = set()
    for q, t in zip(sources, targets):
        if not 0 <= q < len(names) or not 0 <= t < len(names):
            raise ValueError("Edge outside gene universe")
        a, b = names[q], names[t]
        if owners[a] != owners[b]:
            result.add(tuple(sorted((a, b))))
    return result


def component_pairs(names, owners, sources, targets):
    edge_pairs(names, owners, sources, targets)
    graph = coo_matrix((np.ones(len(sources)), (sources, targets)), shape=(len(names), len(names))).tocsr()
    count, labels = connected_components(graph, directed=False)
    groups = [[] for _ in range(count)]
    for name, label in zip(names, labels):
        groups[label].append(name)
    return set(group_pairs(groups, owners))


def transition(before, after, truth):
    gained, lost = after - before, before - after
    return dict(gained_true=len(gained & truth), lost_true=len(lost & truth),
                gained_false=len(gained - truth), lost_false=len(lost - truth))


def analyze_cell(cell, expected_score):
    records = [cell["execution"], cell["native_receipt"], cell["final_partition"], cell["truth"], *cell["inputs"]]
    for item in records:
        check(item)
    execution = json.loads(Path(cell["execution"]["path"]).read_text())
    native = json.loads(Path(cell["native_receipt"]["path"]).read_text())
    artifacts = [execution["private_numeric"], *native["outputs"]]
    for item in artifacts:
        check(item)
    numeric = json.loads(Path(execution["private_numeric"]["path"]).read_text())
    names = numeric["gene_names"]
    owners = dict(zip(names, numeric["gene_to_species"]))
    truth_data = json.loads(Path(cell["truth"]["path"]).read_text())
    truth = {tuple(sorted(pair)) for pair in truth_data["ortholog_pairs"]}
    if len(truth) != len(truth_data["ortholog_pairs"]):
        raise ValueError("Duplicate truth pairs")
    directory = Path(cell["final_partition"]["path"]).parent
    pairs = dict(direct_hits=edge_pairs(names, owners, numeric["hit_queries"], numeric["hit_targets"]))
    for stage in ("rbnh", "multipass"):
        with np.load(directory / (stage + "_edges.npz"), allow_pickle=False) as edges:
            q, t = edges["sources"], edges["targets"]
        pairs[stage + "_edges"] = edge_pairs(names, owners, q, t)
        pairs[stage + "_components"] = component_pairs(names, owners, q, t)
    for stage in ("initial", "multipass", "final"):
        pairs[stage] = set(group_pairs(partition(directory / (stage + ".tsv"), names), owners))
    scores = {stage: score_pairs(pairs[stage], truth_data["ortholog_pairs"], owners) for stage in STAGES}
    if scores["final"] != expected_score:
        raise ValueError("Final stage does not reproduce frozen score")
    transitions = {a + "_to_" + b: transition(pairs[a], pairs[b], truth)
                   for a, b in (("initial", "multipass"), ("multipass", "final"))}
    support = {}
    for label, subset in (("true", pairs["final"] & truth), ("false", pairs["final"] - truth)):
        support[label + "_direct_supported"] = len(subset & pairs["direct_hits"])
        support[label + "_without_direct_hit"] = len(subset - pairs["direct_hits"])
    for item in records + artifacts:
        check(item)
    return dict(condition=cell["condition"], seed=cell["seed"], arm=cell["arm"],
                stages=scores, transitions=transitions, final_direct_support=support,
                evidence=records + artifacts)


def aggregate(rows):
    keys = [(r["condition"], r["seed"], r["arm"]) for r in rows]
    expected = {(c, s, a) for c in CONDITIONS for s in SEEDS for a in ("hmm", "diamond")}
    if len(keys) != len(set(keys)) or set(keys) != expected:
        raise ValueError("Require complete paired reporting panel")
    result = {}
    for condition in (*CONDITIONS, "overall"):
        result[condition] = {}
        for arm in ("hmm", "diamond"):
            selected = [r for r in rows if r["arm"] == arm and (condition == "overall" or r["condition"] == condition)]
            result[condition][arm] = dict(
                stages={stage: {m: mean(r["stages"][stage][m] for r in selected)
                                for m in ("tp", "fp", "fn", "precision", "recall", "f1")} for stage in STAGES},
                transitions={name: {key: mean(r["transitions"][name][key] for r in selected)
                                    for key in selected[0]["transitions"][name]} for name in selected[0]["transitions"]},
                final_direct_support={key: mean(r["final_direct_support"][key] for r in selected)
                                      for key in selected[0]["final_direct_support"]})
    return result


def render(result):
    lines = ["# Matched-Graph Stage Trace", "", "Post-hoc descriptive diagnostic; no causal or confirmatory inference.", "",
             "All values below are equal-dataset mean overlap F1 percentages.",
             "Direct hits, graph edges and component closure are not native ortholog predictions.", "",
             "| Condition | Stage | HMM | DIAMOND |", "|---|---|---:|---:|"]
    for condition, arms in result["summary"].items():
        for stage in STAGES:
            a, b = (100 * arms[arm]["stages"][stage]["f1"] for arm in ("hmm", "diamond"))
            lines.append(f"| {condition} | {stage} | {a:.4f} | {b:.4f} |")
    lines += ["", "## Limits", "", *["- " + s for s in result["limitations"]], ""]
    return "\n".join(lines)


def run(readback_path, scores_path, protocol, output):
    if output.exists():
        raise FileExistsError(output)
    readback_record, scores_record = record(readback_path), record(scores_path)
    if readback_record["sha256"] != READBACK_SHA or scores_record["sha256"] != SCORES_SHA:
        raise ValueError("Unrecognized frozen evidence")
    readback, scores = json.loads(readback_path.read_text()), json.loads(scores_path.read_text())
    check(readback["source"])
    indexed = {(r["condition"], r["seed"], r["arm"]): r for r in scores["records"]}
    rows = []
    for cell in readback["cells"]:
        expected = indexed[cell["condition"], cell["seed"], cell["arm"]]
        if expected["truth"] != cell["truth"] or expected["prediction"] != cell["final_partition"]:
            raise ValueError("Frozen score evidence differs")
        rows.append(analyze_cell(cell, expected["score"]))
    result = dict(status="complete_posthoc_stage_trace", records=rows, summary=aggregate(rows),
                  source=record(__file__), readback=readback_record, scores=scores_record, protocol=record(protocol),
                  dependencies=[record(Path(__file__).with_name(name)) for name in
                                ("simulation_conditions.py", "audit_matched_graph.py", "score_matched_graph.py")],
                  independent_validation=False, inference_rerun=False, scientific_defaults_changed=False,
                  limitations=["Post-hoc, development-exposed diagnostic, not a new confirmatory comparison.",
                               "Overlap scores for search hits, edges and component closure are not ortholog prediction accuracy.",
                               "A pair without direct search evidence can be recovered through indirect graph paths.",
                               "Stage localization does not isolate hit identity, score, ranking or clustering effects causally.",
                               "Five seeds per condition share histories across conditions; no new intervals or tests."])
    output.mkdir(parents=True)
    (output / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    (output / "results.md").write_text(render(result))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--readback", type=Path, required=True)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.readback.resolve(), args.scores.resolve(), args.protocol.resolve(), args.output.absolute())
