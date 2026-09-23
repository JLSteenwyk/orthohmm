"""Compare every NHX line candidate with native reference truth; never choose a tree."""

import argparse
import hashlib
import json
from pathlib import Path

import dendropy
from Bio.Phylo.PhyloXML import Clade, Events, Phylogeny
from benchmark_tools.audit_qfo_swiss_counts import read_raw
from benchmark_tools.compare_swisstree_explicit_events import compare_tree, COUNTS_SHA, MODELS_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen


def candidates(text, members, truth):
    rows = []
    for line_number, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        source = line.strip()
        appended = not source.endswith(";")
        tree = dendropy.Tree.get(data=source + (";" if appended else ""), schema="newick",
            preserve_underscores=True, extract_comment_metadata=True)
        tags = []

        def convert(node):
            flags = [a.value for a in node.annotations if a.name == "D"]
            if len(flags) > 1 or (flags and flags[0] not in ("Y", "N", "T", "F", "?")):
                raise ValueError("Ambiguous or unsupported NHX D tag")
            flag = flags[0] if flags else None
            if node.annotations.get_value("Ev") is not None:
                raise ValueError("NHX Ev semantics require separate review")
            event = Events(duplications=1) if flag in ("Y", "T") else Events(speciations=1) if flag in ("N", "F") else None
            if flag is not None:
                tags.append(dict(value=flag, leaves=sorted(n.taxon.label for n in node.leaf_iter())))
            return Clade(name=node.taxon.label if node.is_leaf() else node.label,
                         events=event, clades=[convert(child) for child in node.child_node_iter()])

        # Transfer the stored hierarchy only; no rooting operation or event inference.
        converted = Phylogeny(root=convert(tree.seed_node), rooted=bool(tree.is_rooted))
        result = compare_tree(converted, members, truth)
        result["rooted_attribute"] = tree.is_rooted
        rows.append(dict(source_line=line_number, source_line_sha256=hashlib.sha256(line.encode()).hexdigest(),
            appended_terminator_for_parse=appended, raw_D_tags=tags, **result))
    if not rows:
        raise ValueError("No NHX candidates")
    comparison_keys = ("pair_counts", "exact_suffix_mapping_candidates", "missing_mapping_candidates",
                       "ambiguous_mapping_candidates", "disagreements")
    return dict(candidate_count=len(rows), candidates=rows, chosen_candidate=None,
        identical_reference_observations=all(all(row[k] == rows[0][k] for k in comparison_keys) for row in rows),
        ancestral_labels_admitted=False)


def compare(root, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    counts_path = root / "benchmark_tools/results/qfo_fastoma_swiss_uncertainty_22098.json"
    models_path = root / "benchmark_tools/results/swisstree_model_acquisition_20260923.json"
    counts = read_frozen(counts_path, COUNTS_SHA)["reconstructed_counts"]
    models = read_frozen(models_path, MODELS_SHA)
    records = [record(counts_path), record(models_path), record(__file__),
        *[record(Path(__file__).with_name(name)) for name in
          ("compare_swisstree_explicit_events.py", "audit_qfo_swiss_counts.py")]]
    reference, members, methods = None, None, []
    for method in counts["methods"]:
        if "raw_file" not in method:
            continue
        item = method["raw_file"]
        check(item)
        records.append(item)
        _, labels, genes = read_raw(Path(item["path"]), counts["families"])
        if reference is None:
            reference, members = labels, genes
        elif labels != reference or genes != members:
            raise ValueError("Native reference truth differs")
        methods.append(method["method"])
    if len(methods) != 7 or len(reference) != counts["reference_relation_count"]:
        raise ValueError("Incomplete native reference panel")
    results = {}
    for identifier, source in models["families"].items():
        if source["tree_format"] != "nhx":
            continue
        if source["status"] != "downloaded_pending_tree_and_mapping_validation":
            raise ValueError("Unacquired NHX source")
        for download in source["downloads"]:
            check(download["file"])
            records.append(download["file"])
        family = source["qfo_family"]
        path = Path(source["downloads"][1]["file"]["path"])
        truth = {(a, b): label for (f, a, b), label in reference.items() if f == family}
        results[family] = dict(identifier=identifier, **candidates(path.read_text(), members[family], truth))
    if set(results) != {"POP", "NOX", "VATB", "SUMF", "APP"}:
        raise ValueError("Changed NHX family inventory")
    for item in records:
        check(item)
    report = dict(status="all_nhx_candidates_compared_not_selected", families=results, records=records,
        native_reference_methods_checked=methods, dendropy_version=dendropy.__version__,
        prediction_statistics_evaluated=False, ancestral_labels_admitted=False, publication_ready=False,
        limitations=["D=Y/T and D=N/F are interpreted conditionally using NHX conventions, not inferred from absent tags.",
            "All physical-line candidates remain retained, even if observations disagree; none is selected by reference agreement.",
            "Stored hierarchy is not independently established biological rooting; parser None/False is preserved.",
            "Unique suffix matches are only candidate aliases, not sequence/taxon identity validation.",
            "Reference relationships share the curated source and are not independent generalization or method-performance evidence."])
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
