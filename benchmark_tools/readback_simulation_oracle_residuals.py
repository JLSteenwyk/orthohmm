"""Independent XML/partition/rule readback of the complete residual pair cohort."""

import argparse
from collections import Counter, defaultdict
from itertools import combinations
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from Bio import Phylo

from benchmark_tools.readback_simulation_gene_tree_oracle import load, record, validate_counts


ORACLE_SHA = "d2299d73238f1a4bb511320add9a62dbb427bc1ea9b4f414e74aefab20e7ad53"
METHOD = "orthohmm_satellite_v2"
CONDITIONS = {"baseline", "divergent", "turnover", "divergent_turnover",
              "missing20", "uneven_taxa", "taxon_count_control"}


def xml_lineages(path):
    root = ET.parse(path).getroot()
    clades = root.findall("./phylogeny/clade")
    if root.tag != "recGeneTree" or len(clades) != 1:
        raise ValueError("XML needs a unique reconciled root")
    paths, descendants, events = {}, {}, {}
    codes = {"speciation": "S", "duplication": "D", "P": "F", "loss": "L"}
    def visit(clade, ancestors):
        name = clade.findtext("name")
        annotations = clade.findall("./eventsRec/*")
        if not name or name in events or len(annotations) != 1 or annotations[0].tag not in codes:
            raise ValueError("Invalid original XML event")
        event = codes[annotations[0].tag]
        events[name] = event
        children = clade.findall("clade")
        if len(children) != (2 if event in {"S", "D"} else 0):
            raise ValueError("Original XML event arity mismatch")
        lineage = (*ancestors, name)
        members = set().union(*(visit(c, lineage) for c in children))
        if event == "F":
            paths[name] = lineage
            members.add(name)
        descendants[name] = members
        return members
    visit(clades[0], ())
    return paths, descendants, events


def owners_from_fasta(items, verify):
    owners = {}
    for item in items:
        path = verify(item.get("absolute_path") or item["path"])
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                gene = line[1:].split()[0]
                if gene in owners:
                    raise ValueError("Duplicate native FASTA header")
                owners[gene] = path.stem
    return owners


def check_candidate(row, expected, paths, descendants, events, owners, species,
                    native_events, active, native_genes, reference):
    genes = set(row["genes"])
    ancestor = row["ancestor"]
    if (genes != native_genes or len(genes) != len(row["genes"])
            or row["status"] != expected["status"] or row["ancestor"] != expected["ancestral_families"][0]
            or row["native_constraint_policy_active"] != active):
        raise ValueError("Candidate/ownership/eligibility differs from original retained inputs")
    labels = {f"F{ancestor}__{g}": g for g in paths}
    if not genes <= labels.keys():
        raise ValueError("Candidate absent from original XML leaves")
    native_clades = {frozenset(f"F{ancestor}__{g}" for g in members) & genes for members in descendants.values()}
    native_clades.discard(frozenset())
    nodes = {n["node_id"]: n for n in row["reconciliation_nodes"]}
    if len(nodes) != len(row["reconciliation_nodes"]):
        raise ValueError("Duplicate reported reconciliation node")
    members = {key: frozenset(n["genes"]) for key, n in nodes.items()}
    if len(set(members.values())) != len(nodes) or set(members.values()) != native_clades:
        raise ValueError("Reported clades differ from independently induced XML clades")
    children = defaultdict(list)
    roots = []
    species_tips = {n.name for n in species.get_terminals()}
    if not {owners[g] for g in genes} <= species_tips:
        raise ValueError("Species tree lacks candidate owners")
    mapped = {key: species.common_ancestor([owners[g] for g in group]) for key, group in members.items()}
    for key, group in members.items():
        supersets = [other for other, part in members.items() if group < part]
        parent = min(supersets, key=lambda other: len(members[other])) if supersets else None
        if nodes[key]["parent_node_id"] != parent:
            raise ValueError("Reported parent differs from induced XML hierarchy")
        children[parent].append(key)
        if parent is None:
            roots.append(key)
    if len(roots) != 1 or members[roots[0]] != genes:
        raise ValueError("Induced XML hierarchy has incorrect root")
    calls = {}
    for key, group in members.items():
        node = nodes[key]
        child_keys = children[key]
        if not child_keys:
            if len(group) != 1 or node["event"] != "leaf" or node["pair_event"] != "leaf":
                raise ValueError("Invalid reported leaf")
            continue
        if len(child_keys) != 2:
            raise ValueError("Induced tree is not binary")
        left, right = [members[c] for c in child_keys]
        overlap = {owners[g] for g in left} & {owners[g] for g in right}
        conflict = any(mapped[c] is mapped[key] for c in child_keys)
        pair_event = "duplication" if overlap else "uncertain" if conflict else "speciation"
        confidence = "high" if not overlap and not conflict or len(overlap) >= 2 else "medium"
        if (node["species_overlap_count"] != len(overlap) or node["mapping_conflict"] != conflict
                or node["pair_event"] != pair_event or node["event_confidence"] != confidence
                or node["event"] != ("duplication" if overlap or conflict else "speciation")
                or node["branch_support"] is not None):
            raise ValueError("Reported calls differ from independent frozen-rule calculation")
        calls[key] = (pair_event, confidence, mapped[key] is species.root and bool(overlap))
    def groups(key):
        descendants_keys = [k for k, group in members.items() if group <= members[key] and k in calls]
        if any(calls[k][2] for k in descendants_keys):
            return [part for child in children[key] for part in groups(child)]
        return [set(members[key])]
    root_groups = groups(roots[0])
    if sorted(map(sorted, root_groups)) != sorted(row["root_groups"]):
        raise ValueError("Reported root groups differ from independently applied species-overlap rule")
    all_pairs = {tuple(sorted(p)) for p in combinations(genes, 2) if owners[p[0]] != owners[p[1]]}
    by_pair_node = {p: min((k for k, group in members.items() if set(p) <= group), key=lambda k: len(members[k]))
                    for p in all_pairs}
    bypass = row["status"] == "unambiguous_bypass"
    raw = all_pairs if bypass else {p for p in all_pairs if calls[by_pair_node[p]][0] != "duplication"}
    high = {p for p in all_pairs if calls[by_pair_node[p]][0] == "speciation"}
    detach, evidence = {}, []
    if active and not bypass:
        for event_index, event in native_events:
            source, target = map(set, (event["source_genes"], event["target_genes"]))
            if not source or not target or source & target or not (source | target) <= genes:
                raise ValueError("Invalid original constraint boundary")
            supporting = next((p for p in sorted(high) if
                p[0] in source and p[1] in target or p[1] in source and p[0] in target), None)
            evidence.append({"event_index": event_index, "supported": supporting is not None,
                             "supporting_pair": list(supporting) if supporting else None,
                             "source_genes": sorted(source), "target_genes": sorted(target)})
            if supporting is None:
                for gene in source:
                    if gene in detach:
                        raise ValueError("Repeated detached source")
                    detach[gene] = event_index
    if evidence != row["constraint_evidence"]:
        raise ValueError("Reported constraint evidence differs from native events/high pairs")
    root_index = {g: i for i, group in enumerate(root_groups) for g in group}
    final_index = {g: (root_index[g], detach.get(g, -1)) if active and not bypass else 0 for g in genes}
    predicted = {p for p in raw if final_index[p[0]] == final_index[p[1]]}
    counts = {"tp": len(predicted & reference), "fp": len(predicted - reference), "fn": len(reference - predicted)}
    if counts != row["counts"] or counts != {k: expected["arms"]["generating_root"][k] for k in counts}:
        raise ValueError("Pair predictions do not reproduce retained oracle counts")
    seen, classes = set(), Counter()
    for pair_row in row["pairs"]:
        pair = (pair_row["gene_a"], pair_row["gene_b"])
        if pair not in all_pairs or pair in seen or pair[0] >= pair[1]:
            raise ValueError("Duplicate, missing or noncanonical pair row")
        seen.add(pair)
        left, right = [labels[g] for g in pair]
        lca = next(a for a, b in reversed(list(zip(paths[left], paths[right]))) if a == b)
        key = by_pair_node[pair]
        truth = events[lca] == "S"
        if truth != (pair in reference) or pair_row["history_node"] != lca or pair_row["history_event"] != events[lca]:
            raise ValueError("Pair ancestor/truth differs from original XML")
        if pair_row["pair_node"] != nodes[key]:
            raise ValueError("Pair uses an incorrect induced reconciliation ancestor")
        # Original XML child clades may contain lost branches; use lineage paths
        # to recover the two surviving branches below this actual ancestor.
        branches = defaultdict(set)
        for leaf, lineage in paths.items():
            if lca in lineage and lineage[-1] != lca:
                branches[lineage[lineage.index(lca) + 1]].add(f"F{ancestor}__{leaf}")
        parts = list(branches.values())
        if len(parts) != 2:
            raise ValueError("Pair ancestor must have two surviving XML branches")
        full_overlap = sorted({owners[g] for g in parts[0]} & {owners[g] for g in parts[1]})
        kept_overlap = sorted({owners[g] for g in parts[0] & genes} & {owners[g] for g in parts[1] & genes})
        flags = {"truth": truth, "predicted": pair in predicted, "raw_predicted": pair in raw,
                 "same_root_group": root_index[pair[0]] == root_index[pair[1]],
                 "same_final_group": final_index[pair[0]] == final_index[pair[1]],
                 "root_filter_active": active and not bypass,
                 "parent_species_overlap": full_overlap, "candidate_species_overlap": kept_overlap,
                 "detachment_events": {g: detach[g] for g in pair if g in detach}}
        if any(pair_row[k] != v for k, v in flags.items()):
            raise ValueError("Pair stage flags or overlap differ from independent readback")
        category = None
        if truth and pair not in predicted:
            category = ("pair_rule_exclusion" if pair not in raw else "root_partition_filter"
                        if root_index[pair[0]] != root_index[pair[1]] else "unsupported_satellite_constraint")
        elif not truth and pair in predicted:
            category = "single_copy_bypass_on_true_duplication" if bypass else "true_duplication_without_retained_species_overlap"
        if pair_row["error_class"] != category:
            raise ValueError("Pair mechanism class differs from independently checked flags")
        if category:
            classes[category] += 1
    if seen != all_pairs or dict(classes) != row["error_classes"]:
        raise ValueError("Incomplete pair cohort or incorrect class projection")
    return {"cell": row["cell"], "family": row["family"], "status": row["status"],
            "pair_rows_verified": len(seen), "counts": counts, "error_classes": dict(classes)}


def run(repo, path, digest):
    report, report_ref = load(path, digest)
    original, _ = load(repo / "benchmark_tools/results/simulation_gene_tree_oracle_20261004.json", ORACLE_SHA)
    expected_cells = {(c, s) for c in CONDITIONS for s in range(20261101, 20261111)}
    if len(original["cells"]) != 70 or {(c["condition"], c["seed"]) for c in original["cells"]} != expected_cells:
        raise ValueError("Incomplete original 70-cell cohort")
    selected = {}
    for cell in original["cells"]:
        if cell["status"] != "baseline_reproduced_oracle_scored":
            raise ValueError("Original cell not complete")
        for candidate in cell["candidates"]:
            score = candidate["arms"]["generating_root"]
            validate_counts(score)
            if score["fp"] or score["fn"]:
                key = (cell["label"], candidate["family"])
                if key in selected:
                    raise ValueError("Duplicate original residual candidate")
                selected[key] = candidate
    rows = report["candidates"]
    if len(rows) != len(selected) or {(r["cell"], r["family"]) for r in rows} != set(selected):
        raise ValueError("Trace does not retain the complete original error cohort")
    inputs = {}
    for item in report["inputs"]:
        if item["path"] in inputs or record(item["path"]) != item:
            raise ValueError("Duplicate or changed selected input")
        inputs[item["path"]] = item
    def verify(path):
        path = Path(path).resolve()
        if str(path) not in inputs:
            raise ValueError("Original readback file not bound in trace inputs")
        return path
    prepared = json.loads(verify(repo / "benchmark_tools/results/simulation_tree_controls_prepared_20260917.json").read_text())
    admission = json.loads(verify(repo / "benchmarks/results/simulation_tree_panel_admission_v1/results.json").read_text())
    native_rows = {r["label"]: r for r in admission["records"] if r["method"] == METHOD and r["variant"] == "generating"}
    datasets = {f"{d['condition']}_{d['seed']}": d for d in prepared["datasets"]}
    checked = []
    for row in rows:
        dataset = datasets[row["cell"]]
        parent = datasets[dataset["parent"]]
        owners = owners_from_fasta(parent["input_evidence"]["inputs"], verify)
        native = native_rows[row["cell"]]
        directory = Path(native["native_report"]["path"]).parent / METHOD
        partition = verify(directory / "orthohmm_working_res/phylogeny_candidate_superfamilies.txt").read_text().splitlines()
        native_genes = set(partition[int(row["family"].removeprefix("Family"))].split())
        reference_doc = json.loads(verify(dataset["input_evidence"]["truth"]["path"]).read_text())
        reference = {tuple(sorted(p)) for p in reference_doc["ortholog_pairs"] if set(p) <= native_genes}
        xml_path = Path(parent["generating_tree"]["path"]).parent.parent / f"G/Gene_trees/{row['ancestor']}_rec.xml"
        paths, descendants, events = xml_lineages(verify(xml_path))
        species = Phylo.read(verify(directory / "orthohmm_phylogeny/species_tree.rooted.nwk"), "newick")
        manifest = json.loads(verify(directory / "orthohmm_phylogeny/provenance_manifest.json").read_text())
        active = manifest["membership_reconciliation"] is not None
        constraint_path = directory / "orthohmm_working_res/phylogeny_candidate_merges.json"
        original_events = json.loads(verify(constraint_path).read_text()) if str(constraint_path) in inputs else []
        local_events = [(i, e) for i, e in enumerate(original_events)
                        if set(e["source_genes"]) & native_genes or set(e["target_genes"]) & native_genes]
        checked.append(check_candidate(row, selected[(row["cell"], row["family"])], paths, descendants,
            events, owners, species, local_events, active, native_genes, reference))
    summary = {"candidates": len(checked), "status_counts": dict(Counter(r["status"] for r in checked)),
               "pair_rows": sum(r["pair_rows_verified"] for r in checked),
               "counts": {k: sum(r["counts"][k] for r in checked) for k in ("tp", "fp", "fn")},
               "error_classes": dict(sum((Counter(r["error_classes"]) for r in checked), Counter()))}
    if (summary != report["summary"] or report["screened_cells"] != 70
            or report["screened_candidates"] != sum(len(c["candidates"]) for c in original["cells"])):
        raise ValueError("Aggregate or screening counts differ from independently checked cohort")
    return {"status": "independent_xml_rules_and_stage_readback_passed", "source": record(__file__),
            "report": report_ref, "input_identities_rechecked": len(inputs), "candidates": checked,
            "summary": summary, "primary_worker_imported": False,
            "limitations": ["Selected post hoc simulation diagnostic, not independent biological validation.",
                "Independently checks XML clades, event ancestors, frozen-rule arithmetic and constraint flags; not native search scores."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = run(args.repo.resolve(), args.report.resolve(), args.report_sha256)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        handle.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result["summary"], sort_keys=True))
