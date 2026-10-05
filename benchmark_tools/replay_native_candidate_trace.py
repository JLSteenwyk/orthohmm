"""Reconstruct retained partitions from accepted merges, not candidate scores."""

import argparse
from collections import defaultdict
import json
from pathlib import Path
import platform
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.diagnose_candidate_trace_variation import key, scalar, trace_comparison, validate_trace
from benchmark_tools.link_factorial_scaling_resources import compare, partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def replay(seed_groups, trace):
    """Apply each round's accepted unions to unchanged round-start memberships."""
    validate_trace(trace)
    current = set(seed_groups)
    if not current or any(not g for g in current):
        raise ValueError("Empty seed partition")
    universe = set().union(*current)
    if sum(map(len, current)) != len(universe):
        raise ValueError("Overlapping seed partition")
    rounds = defaultdict(list)
    for row in trace:
        rounds[row["iteration"]].append(row)
    for iteration in sorted(rounds):
        parents = {g: g for g in current}
        labels = {}

        def find(group):
            while parents[group] != group:
                parents[group] = parents[parents[group]]
                group = parents[group]
            return group

        for row in rounds[iteration]:
            endpoints = []
            for role in ("source", "target"):
                group = frozenset(row[role + "_genes"])
                if group not in current:
                    raise ValueError("Trace endpoint is not a round-start component")
                label = row[role + "_cluster"]
                if label in labels and labels[label] != group:
                    raise ValueError("Conflicting cluster label within round")
                labels[label] = group
                endpoints.append(find(group))
            if endpoints[0] == endpoints[1]:
                raise ValueError("Redundant accepted union")
            parents[endpoints[0]] = endpoints[1]
        components = defaultdict(set)
        for group in current:
            components[find(group)].update(group)
        current = {frozenset(g) for g in components.values()}
    if len(current) != len(seed_groups) - len(trace):
        raise ValueError("Accepted merge count does not explain partition size")
    return universe, current


def collect(preparation_ref, score_ref):
    inputs = {}

    def load(ref):
        check(ref)
        inputs[ref["path"]] = ref
        return Path(ref["path"])

    prep = json.loads(load(preparation_ref).read_text())
    score = json.loads(load(score_ref).read_text())
    if (score.get("schema") != "native_factorial_orthobench_score_v1"
            or score.get("status") != "terminal_native_orthobench_scored"
            or score.get("cell") not in ("p0_c1_r0", "p1_c1_r0")
            or score.get("prediction_format") != "space_separated_groups"
            or score.get("native_outputs_validated") is not True):
        raise ValueError("Require a retained scored candidate-only OrthoBench cell")
    arm = prep["candidate_arms"][score["cell"][:5]]
    if score["original_prediction"] != arm["candidate_partition"]:
        raise ValueError("Original scored partition differs from preparation")
    seed = partition(load(arm["seed_partition"]), "space_separated_groups")
    old = partition(load(arm["candidate_partition"]), "space_separated_groups")
    new = partition(load(score["prediction"]), "space_separated_groups")
    old_trace = json.loads(load(arm["membership_constraints"]).read_text())
    trace_path = str(Path(score["prediction"]["path"]).parent / "phylogeny_candidate_merges.json")
    refs = [ref for ref in score["evidence"] if ref["path"] == trace_path]
    if len(refs) != 1:
        raise ValueError("Native trace is not uniquely bound in retained score evidence")
    new_trace = json.loads(load(refs[0]).read_text())
    for expected, trace in ((old, old_trace), (new, new_trace)):
        if replay(seed[1], trace) != expected:
            raise ValueError("Accepted trace does not reconstruct complete saved partition")
    comparison = compare(old, new, set())
    comparison["native_only_groups"] = comparison.pop("scaling_only_groups")
    # Reference exclusions are inherited from the scored record, not rescored here.
    inherited = score["canonical_partition_comparison"]
    for field in ("reference_genes_in_changed_groups", "reference_touching_groups_equal"):
        comparison[field] = inherited[field]
    if comparison != inherited:
        raise ValueError("Reconstructed partition differences disagree with retained score")
    old_rows, new_rows = validate_trace(old_trace), validate_trace(new_trace)

    def difference(left, right):
        return [{name: scalar(value) for name, value in left[k].items()}
                for k in sorted(left.keys() - right.keys())]

    helpers = [Path(__file__), Path(__file__).with_name("diagnose_candidate_trace_variation.py"),
               Path(__file__).with_name("link_factorial_scaling_resources.py"),
               Path(__file__).with_name("prepare_ob_candidate_neighborhood.py"),
               Path(__file__).with_name("score_ygob_groups.py")]
    for helper in helpers:
        load(record(helper))
    for ref in list(inputs.values()):
        check(ref)
    return {"schema": "native_candidate_accepted_trace_replay_v1", "job_id": score["job_id"],
            "cell": score["cell"], "seed_groups": len(seed[1]), "genes": len(seed[0]),
            "full_original_partition_reconstructed": True, "full_native_partition_reconstructed": True,
            "partition_comparison": comparison,
            "trace_comparison": trace_comparison(old_trace, new_trace,
                arm["expansion"]["parameters"]["max_satellites_per_anchor"]),
            "original_only_accepted_merges": difference(old_rows, new_rows),
            "native_only_accepted_merges": difference(new_rows, old_rows),
            "inputs": sorted(inputs.values(), key=lambda ref: ref["path"]),
            "diagnostic_python": platform.python_version(), "native_inference_reexecuted": False,
            "accuracy_rescored": False, "frozen_method_modified": False,
            "limitations": [
                "Replay of accepted unions reconstructs whole partitions; candidate eligibility, rejected alternatives and score computation are not replayed.",
                "The retained original seed is sufficient for both traces. This does not independently establish every fresh seed's numeric order or the provenance of score-bit differences.",
                "Attachment caps and near ties describe observed decisions, not a demonstrated causal numerical mechanism or a reason to change the frozen method.",
                "Reference exclusions are inherited from the pinned separate score, not new independent accuracy evidence.",
                "Shared-host timing distortion remains unknown and potentially tool-dependent; this diagnostic does not attribute membership variation to contention."]}


def render(report):
    c = report["trace_comparison"]
    p = report["partition_comparison"]
    lines = ["# Accepted Candidate Trace Replay", "",
             f"Job {report['job_id']}, {report['cell']}; read-only retained-data diagnostic.", "",
             f"Both complete partitions reconstruct from {report['seed_groups']:,} original seeds and accepted unions.",
             f"Each has {p['groups']:,} groups covering {p['genes']:,} genes; {p['genes_in_changed_groups']} genes occur in changed groups.",
             f"Common accepted merges: {c['common_semantic_merges']:,} / {c['original_merges']:,}.", "",
             "| Round | Anchor gene | Old / native attachments | Both at cap | Changed-selection support spread |",
             "| ---: | --- | --- | --- | ---: |"]
    for anchor in c["common_anchors_with_changed_sources"]:
        name = next((g for g in anchor["target_genes"] if g.startswith("WBGene")), anchor["target_genes"][0])
        lines.append(f"| {anchor['iteration']} | {name} | {anchor['original_attachments']} / {anchor['native_attachments']} | "
                     f"{anchor['both_at_attachment_cap']} | {anchor['changed_selection_support_spread']:.17g} |")
    lines.extend(["", "Complete changed memberships, selections, numeric differences and direct pins are in diagnostic.json.", "",
                  *["- " + limitation for limitation in report["limitations"]]])
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("preparation", "score", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("preparation-sha256", "score-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    refs = [record(args.preparation), record(args.score)]
    if [r["sha256"] for r in refs] != [args.preparation_sha256, args.score_sha256]:
        raise ValueError("Requested input checksum differs")
    if args.output.exists():
        raise FileExistsError(args.output)
    report = collect(*refs)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "diagnostic.json").write_text(json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n")
    (args.output / "diagnostic.md").write_text(render(report))
    print(json.dumps(record(args.output / "diagnostic.json"), sort_keys=True))
