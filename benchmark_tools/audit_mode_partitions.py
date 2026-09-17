"""Check downstream group partitions in addition to admitted native pair equivalence."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_simulation_generation import verify_file
from benchmark_tools.run_simulation_tree_mode_control import METHOD_SHA
from benchmark_tools.score_ygob_groups import read_predictions, membership
from benchmark_tools.simulation_method_outputs import unique_path

ADMISSION_SHA = "5f3e54e6df5c4e026d03a6bea3245da939e355a0ffdf9a50509cdabf6898e82d"


def canonical(groups, universe):
    index = membership(groups)
    if not set(index) <= universe:
        raise ValueError("Partition has unknown input genes")
    return {tuple(sorted(genes)) for genes in groups.values()}, universe - set(index)


def audit(root, output):
    if output.exists():
        raise FileExistsError(output)
    path = root / "benchmarks/results/simulation_mode_panel_admission_v1/results.json"
    admission = read_frozen(path, ADMISSION_SHA)
    manifest = read_frozen(root / "benchmark_tools/results/publication_variable_native_methods_20260916.json", METHOD_SHA)
    datasets = {d["label"]: d for d in manifest["datasets"]}
    checks = []
    for row in admission["records"]:
        if row["status"] == "unavailable":
            continue
        if row["status"] != "equivalent":
            raise ValueError("Expected admitted equivalent available control")
        method = row["method"]
        if row.get("reused_pilot"):
            check(row["admission"])
            link = json.loads(Path(row["admission"]["path"]).read_text())["pilot_report"]
        else:
            link = row["native_report"]
        check(link)
        native = json.loads(Path(link["path"]).read_text())
        before = Path(datasets[row["label"]]["methods"][method]["output"])
        after = Path(native["provenance"]["configured"]["methods"][method]["output"])
        old_link = native["provenance"]["baseline"][method]["execution_evidence"]
        verify_file(Path(old_link["absolute_path"]), old_link)
        check(native["execution"])
        old_status = json.loads(Path(old_link["absolute_path"]).read_text())
        new_status = json.loads(Path(native["execution"]["path"]).read_text())
        universe = set()
        for item in old_status["verified_inputs"]["inputs"]:
            fasta = Path(item["absolute_path"])
            verify_file(fasta, item)
            for seq in SeqIO.parse(fasta, "fasta"):
                if seq.id in universe:
                    raise ValueError("Duplicate input gene")
                universe.add(seq.id)
        specifications = [("orthogroups", "named_groups", "orthohmm_orthogroups.txt"),
                          ("root_hogs", "root_hogs", "orthohmm_phylogeny/orthohmm_root_hogs.tsv")]
        if method == "orthofinder_full":
            specifications = [("orthogroups", "named_groups", "**/Orthogroups/Orthogroups.txt")]
        for label, format_name, pattern in specifications:
            files = [unique_path(directory, pattern) for directory in (before, after)]
            for file, status in zip(files, (old_status, new_status)):
                matches = [r for r in status["methods"][method]["outputs"] if r["absolute_path"] == str(file)]
                if len(matches) != 1:
                    raise ValueError("Partition absent from original execution inventory")
                verify_file(file, matches[0])
            a, b = [canonical(read_predictions(file, format_name), universe) for file in files]
            checks.append({"dataset": row["label"], "method": method, "partition": label, "input_genes": len(universe),
                           "before": record(files[0]), "after": record(files[1]), "identical": a == b,
                           "groups": [len(a[0]), len(b[0])], "unrepresented_genes": [len(a[1]), len(b[1])],
                           "removed_groups": [list(g) for g in sorted(a[0] - b[0])],
                           "added_groups": [list(g) for g in sorted(b[0] - a[0])]})
    if len(checks) != 199:
        raise ValueError("Incomplete downstream partition inventory")
    read_frozen(path, ADMISSION_SHA)
    result = {"status": "downstream_partitions_checked", "accuracy_evaluated": False, "source": record(__file__),
              "mode_admission": record(path), "checks": checks, "all_identical": all(r["identical"] for r in checks),
              "scope": "Available unchanged-tree controls only; unrepresented gene sets compared explicitly without imputation"}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.root.resolve(), args.output.resolve())
