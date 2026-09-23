"""Describe downloaded model formats and exact-label coverage, not ancestral truth."""

import argparse
import json
from pathlib import Path

from Bio import Phylo
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen

INVENTORY_SHA = "912363276fa5456f11aeff9f398fce71aec8d377c8c8eda7b144105945e554f1"


def inspect_tree(path, format_name, members):
    if format_name not in {"phyloxml", "nhx"}:
        raise ValueError("Unsupported model format")
    tree = Phylo.read(path, "phyloxml" if format_name == "phyloxml" else "newick")
    leaves = [node.name for node in tree.get_terminals()]
    if not leaves or any(not name for name in leaves) or len(set(leaves)) != len(leaves):
        raise ValueError("Missing or duplicate terminal labels")
    candidates = {g: [name for name in leaves if name == g or name.rsplit("_", 1)[-1] == g] for g in members}
    ambiguous = {g: names for g, names in candidates.items() if len(names) > 1}
    matched = {g: names[0] for g, names in candidates.items() if len(names) == 1}
    result = dict(terminal_leaves=len(leaves), leaf_names=leaves,
        exact_label_candidates=matched, ambiguous_label_candidates=ambiguous,
        missing_label_candidates=sorted(g for g, names in candidates.items() if not names),
        benchmark_genes=len(members), mapping_admitted=False,
        rooted_attribute=tree.rooted if format_name == "phyloxml" else None,
        explicit_duplications=None, nodes_with_explicit_duplication_count=None)
    if format_name == "phyloxml":
        events = [node.events.duplications for node in tree.find_clades()
                  if node.events is not None and node.events.duplications is not None]
        if any(type(v) is not int or v < 0 for v in events):
            raise ValueError("Invalid explicit event count")
        result.update(explicit_duplications=sum(events), nodes_with_explicit_duplication_count=len(events))
    return result


def inspect(manifest, inventory, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    inputs = [record(manifest), record(inventory), record(__file__)]
    downloads = json.loads(manifest.read_text())
    families = read_frozen(inventory, INVENTORY_SHA)["family_memberships"]
    if downloads["status"] != "acquisition_finished_not_admitted" or {r["qfo_family"] for r in downloads["families"].values()} != set(families):
        raise ValueError("Incomplete acquisition universe")
    rows = {}
    for identifier, row in downloads["families"].items():
        for download in row["downloads"]:
            check(download["file"])
            inputs.append(download["file"])
        if row["status"] == "unavailable":
            rows[identifier] = dict(status="unavailable", qfo_family=row["qfo_family"], error=row["error"])
            continue
        if row["status"] != "downloaded_pending_tree_and_mapping_validation" or len(row["downloads"]) != 2:
            raise ValueError("Unexpected download contract")
        path = row["downloads"][1]["file"]["path"]
        try:
            rows[identifier] = dict(status="parsed_not_admitted", qfo_family=row["qfo_family"],
                tree_format=row["tree_format"], **inspect_tree(path, row["tree_format"], families[row["qfo_family"]]))
        except Exception as error:
            rows[identifier] = dict(status="parse_failed", qfo_family=row["qfo_family"], error=str(error))
    for item in inputs:
        check(item)
    report = dict(status="model_source_format_and_mapping_feasibility", inputs=inputs, families=rows,
        benchmark_labels_admitted=False, prediction_statistics_evaluated=False, publication_ready=False,
        limitations=["Exact-name/suffix candidates are not validated aliases or topology correspondence.",
            "NHX event semantics are not interpreted; missing is not zero duplication.",
            "XML rooted attributes are retained without rerooting; unlabeled nodes are not called speciations.",
            "No native QfO pair-relation comparison or ancestral-feature calculation has run."])
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("manifest", "inventory", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    inspect(args.manifest, args.inventory, args.output)
