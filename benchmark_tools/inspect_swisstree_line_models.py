"""Diagnose multiple unterminated NHX candidates without selecting a reference."""

import argparse
import hashlib
import json
from pathlib import Path

import dendropy
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def inspect_lines(text):
    candidates = []
    for number, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        normalized = line.strip()
        appended = not normalized.endswith(";")
        tree = dendropy.Tree.get(data=normalized + (";" if appended else ""), schema="newick",
                                preserve_underscores=True, extract_comment_metadata=True)
        names = [node.taxon.label for node in tree.leaf_node_iter()]
        if not names or len(names) != len(set(names)):
            raise ValueError("Empty or repeated candidate leaf names")
        clades, duplications, tags = [], [], []
        for node in tree.preorder_node_iter():
            members = sorted(leaf.taxon.label for leaf in node.leaf_iter())
            if not node.is_leaf():
                clades.append(members)
            flag = node.annotations.get_value("D")
            if flag is not None:
                tags.append(dict(members=members, value=flag))
                if flag == "Y":
                    duplications.append(members)
        candidates.append(dict(source_line=number, source_line_sha256=hashlib.sha256(line.encode()).hexdigest(),
            appended_terminator_for_parse=appended, leaves=sorted(names),
            stored_clades=sorted(clades), explicit_D_Y_clades=sorted(duplications),
            D_tags=tags, rooted_parser_value=tree.is_rooted))
    if not candidates:
        raise ValueError("No candidate tree lines")
    first = candidates[0]
    return dict(candidates=candidates, candidate_count=len(candidates),
        identical_leaf_sets=all(r["leaves"] == first["leaves"] for r in candidates),
        identical_stored_clades=all(r["stored_clades"] == first["stored_clades"] for r in candidates),
        identical_D_Y_clades=all(r["explicit_D_Y_clades"] == first["explicit_D_Y_clades"] for r in candidates),
        chosen_candidate=None, benchmark_labels_admitted=False)


def inspect(directory, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    records, families = [record(__file__)], {}
    for identifier in ("ST001", "ST003", "ST005", "ST018"):
        path = directory / identifier / "modeltree.nhx"
        identity = record(path)
        families[identifier] = inspect_lines(path.read_text())
        check(identity)
        records.append(identity)
    result = dict(status="multiline_model_candidate_diagnostic", families=families, records=records,
        dendropy_version=dendropy.__version__, benchmark_labels_admitted=False,
        prediction_statistics_evaluated=False, publication_ready=False,
        limitations=["Each nonempty physical line is tested as a separate candidate; the source bytes are unchanged.",
            "An explicit semicolon is appended only to the parser input when absent.",
            "Stored clade equality is a serialization comparison, not evidence of biological rooting.",
            "D=Y labels are retained but no source candidate, alias mapping or ancestral interpretation is admitted."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    inspect(args.directory, args.output)
