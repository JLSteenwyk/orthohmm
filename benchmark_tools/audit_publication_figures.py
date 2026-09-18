"""Audit the explicitly retained figure manifests without altering their evidence."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record

PANELS = (
    "figures_accuracy_orthomcl_complete_20260916", "figures_publication_method_20260916",
    "figures_orthobench_factorial_20260916", "figures_ygob_frozen_20260916",
    "figures_ob_sequence_search_20260916", "figures_ob_stratified_errors_20260916",
    "figures_ob_parameter_neighborhood_20260916", "figures_species_tree_robustness_20260916",
    "figures_simulation_variable_native_v2_20260916", "figures_simulation_fixed_native_v2_20260916",
    "figures_simulation_tree_robustness_20260917", "figures_wgd_application_20260917",
    "figures_qfo_swiss_comparators_20260917", "figures/swiss_domain_strata_20260917",
    "qfo_factorial_swiss_figure_20260918",
)


def records(value):
    if isinstance(value, dict):
        if "sha256" in value and "path" in value:
            if set(value) != {"path", "bytes", "sha256"}:
                raise ValueError("Unsupported file-record schema")
            yield value
        else:
            for child in value.values():
                yield from records(child)
    elif isinstance(value, list):
        for child in value:
            yield from records(child)


def inspect_manifest(path, repo, tracked):
    identity = record(path)
    data = json.loads(path.read_text())
    outputs = data["outputs"]
    if not outputs or len({item["path"] for item in outputs}) != len(outputs):
        raise ValueError("Empty or duplicate output inventory")
    if any(Path(item["path"]).resolve().parent != path.resolve().parent for item in outputs):
        raise ValueError("Output is outside its figure directory")
    formats = {Path(item["path"]).suffix for item in outputs}
    if not {".png", ".pdf", ".svg"} <= formats:
        raise ValueError("Missing a required figure format")
    checked = {}
    for expected in records(data):
        name = expected["path"]
        if name in checked:
            if checked[name]["expected"] != expected:
                raise ValueError("Conflicting recorded file identities")
            continue
        file = Path(name)
        relative = str(file.resolve().relative_to(repo)) if file.resolve().is_relative_to(repo) else None
        if not file.is_file():
            actual, state = None, "missing"
        else:
            actual = record(file)
            state = "matches" if actual == expected else "changed"
        checked[name] = {"expected": expected, "actual": actual, "status": state,
                         "repository_relative_path": relative, "tracked_in_main_repository": relative in tracked}
    if record(path) != identity:
        raise ValueError("Manifest changed during inspection")
    for row in checked.values():
        if row["actual"] is not None and record(row["actual"]["path"]) != row["actual"]:
            raise ValueError("File changed during inspection")
    return {"manifest": identity, "status": "all_recorded_bytes_match" if all(
        row["status"] == "matches" for row in checked.values()) else "integrity_failure",
        "output_count": len(outputs), "files": list(checked.values())}


def audit(repo):
    tracked = set(subprocess.check_output(["git", "ls-files", "-z"], cwd=repo, text=True).split("\0"))
    panels = []
    for panel in PANELS:
        path = repo / "benchmark_tools/results" / panel / "manifest.json"
        panels.append({"panel": panel, **inspect_manifest(path, repo, tracked)})
    return {"status": "retained_figure_bytes_verified" if all(p["status"] == "all_recorded_bytes_match" for p in panels)
            else "retained_figure_integrity_failure", "publication_ready": False,
            "source": record(__file__), "panels": panels,
            "total_output_records": sum(p["output_count"] for p in panels),
            "untracked_dependencies": sorted({r["expected"]["path"] for p in panels for r in p["files"]
                                               if not r["tracked_in_main_repository"]}),
            "limitations": ["Explicit retained panels only; historical and defective-runtime figures are not silently promoted.",
                "Byte identity is not scientific, numerical, rendering or statistical validation.",
                "Tracked membership does not establish clean committed bytes, licensing or portable execution.",
                "Manifest source inputs are checked, but their transitive raw-data dependencies are not traversed.",
                "Fixed-length simulation panel remains a failure diagnostic, not an admitted OrthoFinder accuracy comparison.",
                "The QfO factorial figure uses original-release SwissTrees results, not corrected-input reruns.",
                "Pending corrected-QfO factorial and matched scaling figures are not certified by this inventory."]}


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
    raise SystemExit(0 if result["status"] == "retained_figure_bytes_verified" else 1)
