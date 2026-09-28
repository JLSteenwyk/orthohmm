"""Compare preserved tree content across two pinned Selectome archives."""

import argparse
from collections import Counter
import json
from pathlib import Path
import zipfile

import dendropy

from benchmark_tools.inspect_selectome_tf7a import literal_rows, NAME as A_NAME, SHA as A_SHA
from benchmark_tools.inspect_selectome_tf7ab import NAME as AB_NAME, SHA as AB_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def signature(nhx):
    tree = dendropy.Tree.get(data=nhx, schema="newick", preserve_underscores=True,
        suppress_leaf_node_taxa=True, extract_comment_metadata=True)
    labels = [n.label for n in tree.leaf_node_iter()]
    if not labels or any(not s for s in labels) or len(set(labels)) != len(labels):
        raise ValueError("Require unique nonempty case-sensitive leaf labels")
    clades, lengths, events, species, annotations, internal_labels = [], [], [], [], [], []
    descendants = {}
    for node in tree.postorder_node_iter():
        leaves = (node.label,) if node.is_leaf() else tuple(sorted(
            leaf for child in node.child_node_iter() for leaf in descendants[child]))
        descendants[node] = leaves
        clades.append(leaves)
        lengths.append((leaves, node.edge_length))
        events.append((leaves, node.annotations.get_value("D")))
        species.append((leaves, node.annotations.get_value("S")))
        annotations.append((leaves, sorted((a.name, a.value) for a in node.annotations)))
        if not node.is_leaf():
            internal_labels.append((leaves, node.label))
    canonical = lambda rows: sorted(json.dumps(row, sort_keys=True) for row in rows)
    return dict(leaves=sorted(labels), topology=canonical(clades), lengths=canonical(lengths),
        duplication_events=canonical(events), species=canonical(species),
        annotations=canonical(annotations), internal_labels=canonical(internal_labels))


def trees(path, wanted=None):
    found = {}
    with zipfile.ZipFile(path) as z:
        if len(z.infolist()) != 1 or z.namelist()[0] != path.name[:-4]:
            raise ValueError("Unexpected SQL archive member")
        with z.open(z.namelist()[0]) as handle:
            for raw in handle:
                if raw.startswith(b"INSERT INTO `selectome_subtrees` VALUES "):
                    for family, taxon, number, nhx in literal_rows(raw.decode(), "selectome_subtrees"):
                        key = (family, taxon, number)
                        if wanted is None or key in wanted:
                            if key in found:
                                raise ValueError("Duplicate subtree key")
                            found[key] = nhx
    return found


def run(directory):
    refs = [record(directory / name) for name in (A_NAME, AB_NAME)]
    if [r["sha256"] for r in refs] != [A_SHA, AB_SHA]:
        raise ValueError("Pinned archive changed")
    old = trees(directory / A_NAME)
    new = trees(directory / AB_NAME, old)
    if new.keys() != old.keys():
        raise ValueError("Missing shared subtree")
    cases, counts = [], Counter()
    for key in sorted(old):
        a, b = signature(old[key]), signature(new[key])
        differences = [field for field in a if a[field] != b[field]]
        counts.update(differences)
        cases.append(dict(key=key, exact_nhx=old[key] == new[key], differences=differences))
    for ref in refs:
        check(ref)
    return dict(status="shared_selectome_tree_comparison_complete", inputs=refs,
        source=record(__file__), shared_trees=len(cases),
        exact_nhx=sum(c["exact_nhx"] for c in cases), differing_fields=dict(counts),
        cases=cases, benchmark_admitted=False, original_qfo_mapping_recovered=False,
        limitations=["Shared family/taxon/subtree-number keys only; these are derived vertebrate subtrees.",
            "Case-sensitive labels and rooted descendant clades; branch lengths compared exactly.",
            "Duplication events compare D annotations, not an independently inferred biological history.",
            "Other annotations and internal labels are tracked separately; no QfO pair reconstruction."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = run(args.directory)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
