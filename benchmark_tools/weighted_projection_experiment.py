"""Fixed factorial policies and authenticated constituent controls, not tuning."""

import json
import math
from pathlib import Path
import subprocess

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.metric_frontier import ENDPOINTS


PROTOCOL = "hmm_weighted_projection_ablation_20261019.json"
PROTOCOL_REVISION = "7f51f99f7654c9532af20145072897fbe664b703"
LABELS = ("legacy_none", "legacy_henikoff", "legacy_henikoff_cluster", "repair_none",
          "repair_henikoff", "repair_henikoff_cluster")


def specification(source=None):
    root = Path(__file__).resolve().parents[1] if source is None else Path(source)
    return json.loads((root / "benchmark_tools" / PROTOCOL).read_text())


def configuration(label):
    matches = [arm for arm in specification()["arms"] if arm["label"] == label]
    if len(matches) != 1 or label not in LABELS:
        raise ValueError("Unknown weighted projection arm")
    return matches[0]


def stage_signature(graph):
    return [{"label": stage["label"], **{k: stage["output"][k] for k in ("bytes", "sha256")}}
            for stage in graph["stages"]]


def point(row):
    return {"f1": row["scores"][ENDPOINTS[0]], "precision": row["details"]["precision"],
            "recall": row["details"]["recall"]}


def control_references(repo, source):
    """Recover each control from the complete admitted archive, not its score alone."""
    spec = specification(source)
    protocol = source / "benchmark_tools" / PROTOCOL
    if protocol.read_bytes() != subprocess.check_output(
            ["git", "show", PROTOCOL_REVISION + ":benchmark_tools/" + PROTOCOL], cwd=repo):
        raise ValueError("Changed preregistered weighted projection protocol")
    admission = {**spec["individual_admission"], "path": str(repo / spec["individual_admission"]["path"])}
    admission = {k: admission[k] for k in ("path", "bytes", "sha256")}
    check(admission)
    archive = json.loads(Path(admission["path"]).read_text())
    compiler_pin = archive["sources"]["compiler"]["reader"]
    check(compiler_pin)
    compiler = json.loads(Path(compiler_pin["path"]).read_text())
    if (compiler.get("status") != "admitted_individual_hmm_methods_and_projection_transfers_collated"
            or compiler.get("accuracy_result_admitted") is not True
            or compiler.get("defaults_promoted") is not False or compiler.get("goal_completed") is not False
            or tuple(arm["label"] for arm in spec["arms"]) != LABELS):
        raise ValueError("Require the complete admitted individual effects and original factorial design")
    methods = {m["key"]: m for m in compiler["methods"]}
    controls, pins = {}, [record(protocol), admission, compiler_pin]
    for arm in spec["arms"][:4]:
        cells = [c for c in methods[arm["reference_key"]]["cells"] if c["dataset"] == "OrthoBench"]
        if len(cells) != 1 or cells[0]["status"] != "independent_content_readback_verified":
            raise ValueError("Require exactly one admitted constituent control cell")
        cell = cells[0]
        if point(cell) != spec["control_points"][arm["label"]]:
            raise ValueError("Constituent control points differ from preregistration")
        graph_pins = [p for p in cell["checked"] if Path(p["path"]).name == "graph_metrics.json"]
        if len(graph_pins) != 1:
            raise ValueError("Require original authenticated constituent graph")
        graph = json.loads(Path(graph_pins[0]["path"]).read_text())
        native = Path(cell["products"]["orthohmm_root_hogs.tsv"]["path"]).parent
        trees = [p for p in cell["checked"] if Path(p["path"]) == native / "species_tree.rooted.nwk"]
        if len(trees) != 1:
            raise ValueError("Require the original admitted rooted species tree")
        controls[arm["label"]] = {
            "reference_key": arm["reference_key"], "expected": point(cell),
            "candidate": cell["candidate"]["candidate_partition"],
            "constraints": cell["candidate"]["membership_constraints"],
            "golden": {**cell["products"], "species_tree.rooted.nwk": trees[0]},
            "graph_stages": stage_signature(graph), "graph_report": graph_pins[0],
            "baseline_output": str(native.parent), "coverage": cell["coverage"],
        }
        pins += cell["checked"]
    pins = list({p["path"]: p for p in pins}.values())
    for pin in pins:
        check(pin)
    return controls, pins


def validate_activation(graph, label, base_validator):
    arm = configuration(label)
    weighting, projection = arm["profile_sequence_weighting"], arm["profile_projection_policy"]
    result = base_validator(graph, weighting, "weighting", projection_policy=projection)
    records = graph["profile_iterations"][0].get("profile_construction_policies", [])
    joint = arm["role"] == "candidate"
    effective = graph["profile_iterations"][0]["profile_effective_counts"]
    if joint:
        if (len(records) != result["profiles_built"] or not records
                or len({r["cluster_id"] for r in records}) != len(records)
                or any(r["sequence_weighting"] != weighting or r["projection_policy"] != projection
                       or not math.isfinite(r["effective_sequence_count"])
                       or not 0 < r["effective_sequence_count"] <= 200 + 1e-8 for r in records)
                or [{k: r[k] for k in ("cluster_id", "effective_sequence_count")} for r in records] != effective):
            raise ValueError("Missing or inconsistent actual joint construction policy evidence")
    elif records:
        raise ValueError("Constituent control unexpectedly has joint construction records")
    return {**result, "arm": label, "projection_policy": projection,
            "joint_profiles": len(records), "role": arm["role"]}


def validate_control_graph(graph, reference):
    if stage_signature(graph) != reference["graph_stages"]:
        raise ValueError("Constituent graph-stage products differ from their admitted reference")
    return True


def matched_gate(label, observed, controls, gate):
    arm = configuration(label)
    if arm["role"] != "candidate":
        raise ValueError("Constituent controls cannot become advancing candidates")
    matches = {key: gate(controls[key], observed) for key in arm["matched_controls"]}
    return {"eligible": all(g["eligible"] for g in matches.values()), "matched_controls": matches,
            "scope": "All three prespecified development point gates; not significance or independent confirmation"}
