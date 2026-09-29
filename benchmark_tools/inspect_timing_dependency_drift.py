"""Read-only investigation of pyparsing drift in the shared timing runtime."""

import argparse
import importlib.metadata
import json
from pathlib import Path
import sys

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen


def inspect(lookup, sha256, name_diff, output):
    if output.exists():
        raise FileExistsError(output)
    receipt = read_frozen(lookup, sha256)
    differences = json.loads(name_diff.read_text())
    reports, comparisons, imported = [], [], {}
    for method, details in receipt["interpreters"].items():
        ref = details["reports"][-1]
        prior = read_frozen(Path(ref["path"]), ref["sha256"])
        reports.append(ref)
        imported[method] = {k: v for k, v in prior["modules"].items()
                            if k == "pyparsing" or k.startswith("pyparsing.")}
        for old in prior["files"]:
            if "pyparsing" not in Path(old["path"]).parts:
                continue
            current = record(Path(old["path"]))
            comparisons.append(dict(method=method, prior=old, current=current,
                                    same_sha256=old["sha256"] == current["sha256"]))
    distributions = {}
    for name in ("pyparsing", "ProDy", "propka", "dm-tree", "ml-collections"):
        distribution = importlib.metadata.distribution(name)
        metadata = [p for p in distribution.files or [] if str(p).endswith(".dist-info/METADATA")]
        if len(metadata) != 1:
            raise ValueError("Ambiguous installed metadata")
        distributions[name] = dict(version=distribution.version,
                                   metadata=record(Path(distribution.locate_file(metadata[0]))))
    for row in comparisons:
        check(row["current"])
    result = dict(status="shared_timing_runtime_dependency_drift_confirmed",
        sources=[record(lookup), record(name_diff), record(Path(__file__).resolve()), *reports],
        added_entry_count=len(differences["added"]), removed_entry_count=len(differences["removed"]),
        entry_name_differences=differences, historical_imports=imported,
        current_distributions=distributions, file_comparisons=comparisons,
        changed_imported_files=sum(not r["same_sha256"] for r in comparisons),
        python=sys.executable, scientific_execution_authorized=False,
        limitations=["Read-only diagnosis; no dependency restoration, substitution, native fixture or timing run.",
                     "Path-name inventory is not a full content audit; changed pyparsing files are individually hash checked.",
                     "Does not identify who changed the environment or establish an effect on scientific outputs."])
    with output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lookup", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--name-diff", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = inspect(args.lookup.resolve(), args.sha256, args.name_diff.resolve(), args.output.resolve())
    print(json.dumps({k: result[k] for k in ("added_entry_count", "removed_entry_count", "changed_imported_files")}))
