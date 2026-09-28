"""Inventory the checksummed Selectome TF7A dump without executing SQL."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import zipfile

import dendropy

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

NAME = "selectome_03-TF7A__mysql5.0.sql.zip"
SHA = "4047be74d730e7057633b5282abfe999cc0202dd42eb194b8efc69586c2745d1"


def literal_rows(statement, table):
    import sqlglot
    from sqlglot import exp

    expressions = sqlglot.parse(statement, read="mysql")
    if len(expressions) != 1:
        raise ValueError("Require one INSERT statement")
    node = expressions[0]
    if (not isinstance(node, exp.Insert) or not isinstance(node.this, exp.Table)
            or node.this.name != table or not isinstance(node.expression, exp.Values)):
        raise ValueError("Unexpected table or nonliteral INSERT")
    for row in node.expression.expressions:
        if not isinstance(row, exp.Tuple) or any(not isinstance(v, exp.Literal) for v in row.expressions):
            raise ValueError("Only literal tuples are accepted")
        yield [v.to_py() for v in row.expressions]


def inspect(directory):
    import sqlglot

    archive, checksums = directory / NAME, directory / "SHA256.sum"
    refs = [record(archive), record(checksums)]
    published = [line.split() for line in checksums.read_text().splitlines()]
    if refs[0]["sha256"] != SHA or [r for r in published if r[-1] == NAME] != [[SHA, NAME]]:
        raise ValueError("Archive differs from pinned published checksum")
    with zipfile.ZipFile(archive) as z:
        entries = z.infolist()
        if len(entries) != 1 or entries[0].filename != NAME[:-4] or entries[0].file_size != 15917948:
            raise ValueError("Unexpected archive inventory")
        if z.testzip() is not None:
            raise ValueError("ZIP CRC failure")
        payload = z.read(entries[0])
    data = payload.decode("utf-8")
    subtrees, news, tables = [], [], []
    keys, families, species = set(), set(), Counter()
    for line in data.splitlines():
        if line.startswith("CREATE TABLE "):
            tables.append(line.split("`")[1])
        if line.startswith("INSERT INTO `selectome_news` VALUES "):
            news.extend(literal_rows(line, "selectome_news"))
        if not line.startswith("INSERT INTO `selectome_subtrees` VALUES "):
            continue
        for family, taxon, number, nhx in literal_rows(line, "selectome_subtrees"):
            key = (family, taxon, number)
            if key in keys or not isinstance(nhx, str) or "&&NHX" not in nhx:
                raise ValueError("Duplicate subtree or missing NHX content")
            tree = dendropy.Tree.get(data=nhx, schema="newick", preserve_underscores=True,
                                     suppress_leaf_node_taxa=True,
                                     extract_comment_metadata=True)
            leaves = list(tree.leaf_node_iter())
            labels = Counter(leaf.label for leaf in leaves)
            owners = [leaf.annotations.get_value("S") for leaf in leaves]
            if not leaves or any(not s for s in owners):
                raise ValueError("Missing leaf species annotations")
            keys.add(key)
            families.add(family)
            species.update(owners)
            subtrees.append(dict(family=family, taxon=taxon, number=number,
                leaves=len(leaves), species=sorted(set(owners)),
                duplicate_leaf_labels={label: count for label, count in labels.items() if count > 1},
                nhx_sha256=hashlib.sha256(nhx.encode()).hexdigest()))
    if not subtrees or not any("TreeFam-A 7" in row[0] for row in news):
        raise ValueError("Missing tree records or release identification")
    for ref in refs:
        check(ref)
    return dict(status="selectome_tf7a_inventory_complete", inputs=refs,
        source=record(__file__), sqlglot_version=sqlglot.__version__,
        dendropy_version=dendropy.__version__, archive_crc_passed=True,
        sql_member_sha256=hashlib.sha256(payload).hexdigest(), tables=tables,
        news=news, subtrees=subtrees, subtree_count=len(subtrees), family_count=len(families),
        leaf_occurrences_by_species=dict(species),
        mapping_filename_present="treefam2reference" in data,
        original_qfo_inputs_recovered=False, benchmark_admitted=False,
        limitations=["Selectome-derived subtrees, not authenticated complete QfO source trees.",
            "No QfO reference-generation mapping recovered; no family labels assigned to benchmark pairs.",
            "SQL parsed as data only; no database import or downloaded SQL execution."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = inspect(args.directory)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
