"""Compare explicit stored-tree events with native reference labels, never method scores."""

import argparse
from collections import Counter
import json
from pathlib import Path

from Bio import Phylo
from benchmark_tools.audit_qfo_swiss_counts import read_raw
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

COUNTS_SHA = "121cc8adf3cd63878f19006ec5500a13f879d042eccd565e1ae5f5a863f43fb1"
MODELS_SHA = "51baa95224887eb25e7c4520d8c8cf0277100c32971dfee1a242e58d32473e61"


def event_label(node):
    event = node.events
    if event is None:
        return "unknown"
    counts = (event.duplications, event.speciations)
    if any(v is not None and (type(v) is not int or v < 0) for v in counts):
        raise ValueError("Invalid explicit event count")
    duplication = bool(event.duplications)
    speciation = bool(event.speciations)
    if duplication and speciation:
        return "ambiguous"
    return "duplication" if duplication else "speciation" if speciation else "unknown"


def compare_tree(tree, members, truth):
    leaves = tree.get_terminals()
    names = [leaf.name for leaf in leaves]
    if not names or any(not name for name in names) or len(set(names)) != len(names):
        raise ValueError("Missing or duplicate tree leaf labels")
    candidates = {gene: [leaf for leaf in leaves if leaf.name == gene or leaf.name.rsplit("_", 1)[-1] == gene]
                  for gene in members}
    mapping = {gene: nodes[0] for gene, nodes in candidates.items() if len(nodes) == 1}
    if len({id(node) for node in mapping.values()}) != len(mapping):
        raise ValueError("Different accessions map to the same tree leaf")
    counts = Counter()
    disagreements = []
    for (a, b), ortholog in sorted(truth.items()):
        if a not in members or b not in members or a == b or type(ortholog) is not bool:
            raise ValueError("Invalid reference pair")
        if a not in mapping or b not in mapping:
            counts["unmapped_pairs"] += 1
            continue
        node = tree.common_ancestor(mapping[a], mapping[b])
        event = event_label(node)
        counts[event + ("_reference_ortholog" if ortholog else "_reference_nonortholog")] += 1
        if (event == "duplication" and ortholog) or (event == "speciation" and not ortholog):
            disagreements.append(dict(a=a, b=b, stored_lca_event=event, reference_ortholog=ortholog,
                                      stored_clade_leaves=sorted(n.name for n in node.get_terminals())))
    return dict(rooted_attribute=tree.rooted, reference_pairs=len(truth), pair_counts=dict(sorted(counts.items())),
        exact_suffix_mapping_candidates={g: node.name for g, node in sorted(mapping.items())},
        missing_mapping_candidates=sorted(g for g, nodes in candidates.items() if not nodes),
        ambiguous_mapping_candidates={g: [n.name for n in nodes] for g, nodes in candidates.items() if len(nodes) > 1},
        disagreements=disagreements, mapping_admitted=False, ancestral_labels_admitted=False)


def compare(root, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    count_path = root / "benchmark_tools/results/qfo_fastoma_swiss_uncertainty_22098.json"
    model_path = root / "benchmark_tools/results/swisstree_model_acquisition_20260923.json"
    counts = read_frozen(count_path, COUNTS_SHA)["reconstructed_counts"]
    models = read_frozen(model_path, MODELS_SHA)
    records = [record(__file__), record(Path(__file__).with_name("audit_qfo_swiss_counts.py")),
               record(count_path), record(model_path)]
    reference, members, checked_methods = None, None, []
    for method in counts["methods"]:
        if "raw_file" not in method:
            continue
        raw = method["raw_file"]
        check(raw)
        records.append(raw)
        _, labels, membership = read_raw(Path(raw["path"]), counts["families"])
        if reference is None:
            reference, members = labels, membership
        elif labels != reference or membership != members:
            raise ValueError("Reference labels differ between native assessments")
        checked_methods.append(method["method"])
    if len(checked_methods) != 7 or len(reference) != counts["reference_relation_count"]:
        raise ValueError("Incomplete native reference comparison")
    if {row["qfo_family"] for row in models["families"].values()} != set(members):
        raise ValueError("Tree and reference family inventories differ")
    results = {}
    for identifier, source in models["families"].items():
        family = source["qfo_family"]
        if source["status"] != "downloaded_pending_tree_and_mapping_validation":
            raise ValueError("Incomplete model acquisition")
        for item in source["downloads"]:
            check(item["file"])
            records.append(item["file"])
        if source["tree_format"] != "phyloxml":
            results[family] = dict(status="nhx_semantics_not_compared", identifier=identifier)
            continue
        path = Path(source["downloads"][1]["file"]["path"])
        tree = Phylo.read(path, "phyloxml")
        truth = {(a, b): label for (f, a, b), label in reference.items() if f == family}
        results[family] = dict(status="stored_xml_event_diagnostic_only", identifier=identifier,
                              **compare_tree(tree, members[family], truth))
    for item in records:
        check(item)
    report = dict(status="explicit_stored_xml_events_compared_not_admitted", records=records,
        native_reference_methods_checked=checked_methods, families=results,
        prediction_statistics_evaluated=False, ancestral_labels_admitted=False, publication_ready=False,
        limitations=["Only reference truth (TP/FN versus FP/TN) is retained; no method-performance comparison or feature bins.",
            "Stored orientation is used diagnostically; rooted=false is not silently changed or treated as rooted biological history.",
            "Absent event annotations are unknown, never inferred speciations; NHX candidates are not selected.",
            "Exact-name/suffix matches are mapping candidates, not independently validated accession/taxon/sequence aliases.",
            "Current curated trees may differ in version and sampling from QfO 2020; agreement is not independent test evidence."])
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    compare(args.root.resolve(), args.output.absolute())
