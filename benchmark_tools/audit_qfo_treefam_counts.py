"""Validate pooled TreeFam-A counts without inventing independent families."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_swiss_counts import read_raw, statistics, SCORER_SHA, IMAGE_SHA
from benchmark_tools.report_qfo_recovered_stages import ADMISSION_SHA, STAGES
from benchmark_tools.snapshot_orthohmm_input_order import record

REFERENCE_SHA = "6f419f96886e5cf14ac889a1437cdb2655140cb8fd281f6bb5f8e0baa69d23b3"


def verify_counts(counts, native):
    scores = statistics(counts)
    for metric, axis in (("TPR", "metric_x"), ("PPV", "metric_y")):
        if not math.isclose(scores[metric], native["native_participant"][axis], rel_tol=0, abs_tol=5e-8):
            raise ValueError("Pooled counts do not reproduce native " + metric)
    if not math.isclose(scores["F1"], native["score"], rel_tol=0, abs_tol=5e-8):
        raise ValueError("Harmonic mean differs from retained endpoint")
    return scores


def audit(repo):
    checked = []
    def check(path, digest=None):
        item = record(path)
        if digest is not None and item["sha256"] != digest:
            raise ValueError("Changed frozen input: " + str(path))
        checked.append(item)
        return Path(path)
    results = repo / "benchmark_tools/results"
    swiss = json.loads((results / "qfo_swiss_counts_20260917.json").read_text())
    admission_path = check(Path(swiss["admission"]["path"]), ADMISSION_SHA)
    admission = json.loads(admission_path.read_text())
    if admission["status"] != "four_stage_assessments_checked" or [r["stage"] for r in admission["records"]] != list(STAGES):
        raise ValueError("Changed stage admission")
    base = repo / "qfo_benchmark"
    image = check(base / "scoring/container_cache/qfobenchmark-darwin-2022.1.img", IMAGE_SHA)
    reference = check(base / "benchmark-webservice/reference_data/2020/ReconciledTrees_TreeFam-A.drw", REFERENCE_SHA)
    scorer = check(base / "benchmark-webservice/RefPhyloTest.drw", SCORER_SHA)
    generator = check(base / "benchmark-webservice/generateData/AddReconciledTree.drw")
    script = check(Path(__file__).with_name("audit_qfo_treefam_reference.drw"))
    if "'" in str(reference):
        raise ValueError("Unsupported Darwin path quoting")
    command = ["singularity", "exec", str(image), "darwin", "-E"]
    result = subprocess.run(command, input=f"reference := '{reference}':\n" + script.read_text(),
                            text=True, capture_output=True, check=True, timeout=60)
    inventory = re.findall(r"^TREEFAM_REFERENCE\t([^\t\n]+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)$", result.stdout, re.MULTILINE)
    if len(inventory) != 1 or inventory[0][0] != "TreeFamA" or int(inventory[0][3]) != 0 or int(inventory[0][5]) != 0:
        raise ValueError("Pooled one-direction reference assumption failed")
    name, n_members, n_relations, _, n_incident, _ = inventory[0]
    container_source = subprocess.run(["singularity", "exec", str(image), "cat", "/benchmark/RefPhyloTest.drw"],
                                      capture_output=True, check=True, timeout=30).stdout
    if hashlib.sha256(container_source).hexdigest() != SCORER_SHA:
        raise ValueError("Container scorer differs from inspected source")
    stages, expected_truth, expected_members = [], None, None
    for index, row in enumerate(admission["records"]):
        if row["index"] != index or row["status"] != "admitted":
            raise ValueError("Unadmitted stage")
        paths = list((base / f"scoring/checked_v2_{index}/results/TreeFam-A").glob("*raw.txt.gz"))
        if len(paths) != 1:
            raise ValueError("Missing or ambiguous raw TreeFam evidence")
        raw = check(paths[0])
        counts, truth, members = read_raw(raw, [name])
        if sum(counts[name].values()) != int(n_relations) or len(members[name]) != int(n_incident):
            raise ValueError("Raw rows do not cover reference inventory")
        if expected_truth is None:
            expected_truth, expected_members = truth, members
        if truth != expected_truth or members != expected_members:
            raise ValueError("Stage reference truth or member universe changed")
        scores = verify_counts(counts[name], row["assessment"]["endpoints"]["TreeFam-A"])
        stages.append({"stage": row["stage"], "raw": record(raw), "counts_without_prior": dict(counts[name]),
                       "native_statistics_with_prior": scores})
    for item in checked:
        if record(item["path"]) != item:
            raise ValueError("Artifact changed during audit")
    return {"status": "pooled_treefam_counts_verified", "source": record(__file__),
            "shared_count_reader": record(Path(__file__).with_name("audit_qfo_swiss_counts.py")),
            "checked_inputs": checked, "reference_cases": 1, "reference_case": name,
            "reference_members": int(n_members), "reference_relations": int(n_relations),
            "relation_incident_members": int(n_incident), "reference_members_without_relations": int(n_members) - int(n_incident),
            "orientation_evidence": {"command": command, "stdout": result.stdout, "stderr": result.stderr},
            "stages": stages, "count_conversion": "raw one-direction count / 2 + 1, one prior for the entire pooled case",
            "family_level_uncertainty": "not_identifiable_from_retained_case_labels",
            "limitations": ["TreeFam-A has one pooled native case, not one case per original family.",
                            "Original source-family mapping must be recovered and validated before family resampling.",
                            "A single-case bootstrap is degenerate; independent-pair resampling ignores dependence.",
                            "Reference-graph components have not been established as original independent families.",
                            "These are count and reference-consistency checks, not significance or superiority evidence.",
                            "Current raw hashes and generator source are retrospective provenance, not execution-time attestations."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.repo.resolve())
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
