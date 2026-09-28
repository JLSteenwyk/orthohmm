"""Stream the pinned TF7AB SQL dump as data, without importing a database."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import zipfile

import dendropy

from benchmark_tools.inspect_selectome_tf7a import literal_rows
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

NAME = "selectome_04-TF7AB__mysql5.0.sql.zip"
SHA = "d936fa7fd1bbaa96904d787c705f54f2aeab4ad936f3ab4cd9ace41f0595d6e0"


def inspect(directory, earlier):
    import sqlglot

    refs = [record(directory / NAME), record(directory / "SHA256.sum"), record(earlier)]
    checksums = [line.split() for line in (directory / "SHA256.sum").read_text().splitlines()]
    if refs[0]["sha256"] != SHA or [r for r in checksums if r[-1] == NAME] != [[SHA, NAME]]:
        raise ValueError("Published checksum mismatch")
    old = json.loads(earlier.read_text())
    if old["inputs"][0]["sha256"] != "4047be74d730e7057633b5282abfe999cc0202dd42eb194b8efc69586c2745d1":
        raise ValueError("Wrong earlier TF7A source")
    tables, news, schema = [], [], []
    indices, identifiers, families = set(), set(), set()
    taxonomy, subtree_taxa, species = Counter(), Counter(), Counter()
    trees, leaves, exact_duplicates = {}, 0, 0
    digest, mapping_marker = hashlib.sha256(), False
    with zipfile.ZipFile(directory / NAME) as z:
        entries = z.infolist()
        if len(entries) != 1 or entries[0].filename != NAME[:-4] or entries[0].file_size != 371301900:
            raise ValueError("Unexpected archive inventory")
        in_schema = False
        with z.open(entries[0]) as handle:
            for raw in handle:
                digest.update(raw)
                line = raw.decode("utf-8")
                mapping_marker |= "treefam2reference" in line
                if line.startswith("CREATE TABLE "):
                    tables.append(line.split("`")[1])
                    in_schema = tables[-1] == "genes"
                if in_schema:
                    schema.append(line.rstrip())
                    if line.startswith(")"):
                        in_schema = False
                if line.startswith("INSERT INTO `selectome_news` VALUES "):
                    news.extend(literal_rows(line, "selectome_news"))
                elif line.startswith("INSERT INTO `genes` VALUES "):
                    for row in literal_rows(line, "genes"):
                        if len(row) != 10 or row[0] in indices or row[1] in identifiers:
                            raise ValueError("Unexpected gene row or duplicate identity")
                        indices.add(row[0])
                        identifiers.add(row[1])
                        taxonomy[row[8]] += 1
                elif line.startswith("INSERT INTO `selectome_subtrees` VALUES "):
                    for family, taxon, number, nhx in literal_rows(line, "selectome_subtrees"):
                        key = (family, taxon, number)
                        if key in trees or "&&NHX" not in nhx:
                            raise ValueError("Duplicate subtree or missing NHX")
                        tree = dendropy.Tree.get(data=nhx, schema="newick", preserve_underscores=True,
                            suppress_leaf_node_taxa=True, extract_comment_metadata=True)
                        nodes = list(tree.leaf_node_iter())
                        owners = [n.annotations.get_value("S") for n in nodes]
                        if not nodes or any(not s for s in owners):
                            raise ValueError("Missing leaf species")
                        leaves += len(nodes)
                        exact_duplicates += int(len({n.label for n in nodes}) != len(nodes))
                        species.update(owners)
                        families.add(family)
                        subtree_taxa[taxon] += 1
                        trees[key] = hashlib.sha256(nhx.encode()).hexdigest()
        # Reading the sole member to EOF also checks its ZIP CRC.
    if not trees or not identifiers or not any("TreeFam7 A &amp; B" in r[0] for r in news):
        raise ValueError("Missing expected TF7AB contents")
    prior = {(r["family"], r["taxon"], r["number"]): r["nhx_sha256"] for r in old["subtrees"]}
    shared = prior.keys() & trees.keys()
    for ref in refs:
        check(ref)
    return dict(status="selectome_tf7ab_inventory_complete", inputs=refs,
        source=record(__file__), parser_helper=record(Path(__file__).with_name("inspect_selectome_tf7a.py")),
        sqlglot_version=sqlglot.__version__, dendropy_version=dendropy.__version__,
        sql_member_sha256=digest.hexdigest(), archive_crc_passed=True, tables=tables,
        gene_schema=schema, gene_count=len(identifiers), genes_by_taxonomy=dict(taxonomy),
        news=news, subtree_count=len(trees), family_count=len(families), subtree_taxa=dict(subtree_taxa),
        leaf_occurrences=leaves, leaf_occurrences_by_species=dict(species),
        trees_with_exact_duplicate_leaf_labels=exact_duplicates,
        tree_inventory_sha256=hashlib.sha256(json.dumps(sorted(trees.items())).encode()).hexdigest(),
        earlier_tree_comparison=dict(earlier=len(prior), shared_keys=len(shared),
            identical_nhx=sum(prior[k] == trees[k] for k in shared), missing_keys=len(prior.keys()-trees.keys())),
        mapping_filename_present=mapping_marker, original_qfo_inputs_recovered=False,
        benchmark_admitted=False, limitations=["Derived Selectome subtrees, not authenticated complete QfO input trees.",
            "Gene identifiers and descriptions are not the original QfO mapping; no benchmark pairs assigned to families.",
            "No SQL executed; exact NHX comparison is conditional on matching family/taxon/subtree-number keys."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("directory", "earlier", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = inspect(args.directory, args.earlier)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
