"""Explicit prospective inventory adoption; never edit the historical runtime."""

import argparse
from copy import deepcopy
import json
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_controlled_fragment_observations import PINS, record
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, execution_environment

PARENT = "publication_variable_native_methods_20260916.json"
SCHEMA = "controlled_fragment_execution_runtime_v1"
CRITICAL = ("numpy", "numba", "llvmlite", "biopython", "scipy", "igraph", "leidenalg",
            "dendropy", "psutil", "texttable", "packaging")
QUERY = "import importlib.metadata as m,json,sys; names=sorted({d.metadata['Name'] for d in m.distributions() if d.metadata['Name']}); print(json.dumps({'python':sys.version,'packages':{n:m.version(n) for n in names}}))"


def package_version(packages, name):
    matches = [version for key, version in packages.items() if re.sub(r"[-_.]+", "-", key).lower() == name]
    if len(matches) != 1 or not isinstance(matches[0], str) or not matches[0]:
        raise ValueError("Missing or ambiguous required distribution: " + name)
    return matches[0]


def adopt_inventory(parent, current, parent_record):
    """Keep scientific dependencies exact; label all other inventory differences."""
    if set(current) != {"orthohmm", "orthofinder"}:
        raise ValueError("Both tool inventories required")
    if current["orthofinder"] != parent["environments"]["orthofinder"]:
        raise ValueError("OrthoFinder inventory changed; this recovery does not authorize it")
    prior = parent["environments"]["orthohmm"]
    now = current["orthohmm"]
    if prior["python"] != now["python"]:
        raise ValueError("OrthoHMM interpreter version changed")
    critical = {}
    for name in CRITICAL:
        old, new = package_version(prior["packages"], name), package_version(now["packages"], name)
        if old != new:
            raise ValueError("Scientific dependency changed: " + name)
        critical[name] = new
    differences = {name: {"historical": prior["packages"].get(name), "prospective": now["packages"].get(name)}
                   for name in sorted(set(prior["packages"]) | set(now["packages"]))
                   if prior["packages"].get(name) != now["packages"].get(name)}
    result = deepcopy(parent)
    result.update(schema=SCHEMA, status="prospective_fragment_runtime", parent_manifest=deepcopy(parent_record),
                  environments=deepcopy(current), inventory_amendment={
                      "scope": "Explicit new OrthoHMM inventory; historical metadata is not rewritten",
                      "scientific_dependencies_unchanged": critical, "differences": differences,
                      "historical_full_inventory_equal": not differences,
                      "historical_output_equivalence_established": False,
                      "limitations": ["Dependency versions and source/native binary checks are not proof of identical historical outputs.",
                                      "Original package-inventory preflight refusal remains a retained failed check, not an inference attempt."]})
    return result


def freeze(root, output):
    if output.exists():
        raise FileExistsError("Existing prospective runtime; do not replace provenance")
    parent_path = root / "benchmark_tools/results" / PARENT
    parent = read_frozen(parent_path, PINS[PARENT])
    interpreters = {"orthohmm": parent["tool_entrypoints"]["orthohmm_python"]["absolute_path"],
                    "orthofinder": str(Path(parent["tool_entrypoints"]["orthofinder"]["absolute_path"]).parent / "python")}
    current = {name: json.loads(subprocess.check_output([python, "-c", QUERY], text=True))
               for name, python in interpreters.items()}
    result = adopt_inventory(parent, current, record(parent_path))
    verify_environment(result)
    _, resolved = execution_environment(result)
    result.update(adoption_source=record(__file__), path_resolved_executables=resolved)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as handle:
        handle.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = freeze(args.root.resolve(), args.output.resolve())
    print(json.dumps({"status": result["status"], "schema": result["schema"],
                      "inventory_changes": len(result["inventory_amendment"]["differences"]),
                      "scientific_dependencies": result["inventory_amendment"]["scientific_dependencies_unchanged"],
                      "historical_output_equivalence_established": False}, sort_keys=True))
