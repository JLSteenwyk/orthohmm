"""Read actual selective-profile receipts before admitting pilot activation."""

import json
import math
from pathlib import Path

from Bio import SeqIO

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from orthohmm.search.selective_hmmer import read_hmmer_tables
from orthohmm.search.selective_profiles import decide_query, query_list_digest


def validate(graph, mode, folder, gene_names, expected_tools=None):
    if mode == "off":
        if graph.get("selective_profile_evidence") or graph["parameters"].get("selective_profile_mode", "off") != "off":
            raise ValueError("Off control contains selective profile intervention")
        return {"selective_mode": "off"}
    root = Path(folder).resolve() / "graph/selective_profiles"
    receipt = root / "selection.json"
    evidence, = graph["selective_profile_evidence"]
    if (json.loads(receipt.read_text()) != evidence or evidence["mode"] != mode
            or evidence["status"] != "selective_profile_edges_built" or evidence["orthology_inferred"] is not False
            or graph["parameters"].get("selective_profile_mode") != mode):
        raise ValueError("Actual selective receipt differs from graph metadata or expected mode")
    checked = [record(receipt)]
    selection = evidence["selection"]
    selected = selection["selected"]
    if (len(selected) > 128 or selection["candidate_digest"] != query_list_digest(selected)
            or selection["orthology_inferred"] is not False):
        raise ValueError("Invalid fixed candidate list/digest")
    partition = Path(folder) / "graph/orthogroups_multipass.txt"
    groups = [line.split() for line in partition.read_text().splitlines() if line.split()]
    genes = [gene for group in groups for gene in group]
    if len(genes) != len(set(genes)) or set(genes) != set(gene_names):
        raise ValueError("Original profile-source partition has incomplete unique ownership")
    checked.append(record(partition))
    owners = {gene: family for family, group in enumerate(groups) for gene in group}
    by_query, by_family, thresholds = {}, {}, {}
    for query in selected:
        q, name = query["query_id"], query["query_gene"]
        if (type(q) is not int or not 0 <= q < len(gene_names) or gene_names[q] != name or q in by_query
                or query["status"] != "selected" or not math.isfinite(query["priority"]) or query["priority"] < 0):
            raise ValueError("Invalid selected query identity/priority")
        by_query[q] = query
        competitors = query["competitors"]
        external = [c for c in competitors if not c["is_own"]]
        if (not 1 <= len(external) <= 4 or not len(external) <= len(competitors) <= 5
                or len({c["profile_id"] for c in competitors}) != len(competitors)
                or query["external_profiles"] != [c["profile_id"] for c in external]
                or not any(c["anchor"] is not None for c in external)):
            raise ValueError("Incomplete competing-profile set or missing positive support")
        if any(not math.isfinite(c["original_hit_score"]) for c in external):
            raise ValueError("Invalid original external score")
        best = min(external, key=lambda c: (-c["original_hit_score"], tuple(sorted(c["members"]))))
        priority = abs(best["original_hit_score"] - best["original_self_threshold"]) / max(1., abs(best["original_self_threshold"]))
        if (query["priority_profile_id"] != best["profile_id"] or query["priority"] != priority
                or (len(external) == 1 and best["original_hit_score"] >= best["original_self_threshold"])):
            raise ValueError("Ambiguity or original boundary-priority rule differs")
        own = [c["profile_id"] for c in competitors if c["is_own"]]
        if (query["own_profile_id"] not in (None, owners[name])
                or own != ([] if query["own_profile_id"] is None else [owners[name]])):
            raise ValueError("Own-family competitor differs from source ownership")
        for c in competitors:
            family = c["profile_id"]
            if (type(family) is not int or not 0 <= family < len(groups)
                    or c["members"] != groups[family] or not 3 <= len(c["members"]) <= 200
                    or c["is_own"] != (family == owners[name])
                    or not math.isfinite(c["original_self_threshold"])):
                raise ValueError("Profile source members or threshold differ from original family")
            if family in thresholds and thresholds[family] != c["original_self_threshold"]:
                raise ValueError("Inconsistent original source threshold")
            thresholds[family] = c["original_self_threshold"]
            by_family.setdefault(family, []).append(name)
            anchor = c["anchor"]
            if anchor is not None and (
                type(anchor["gene_id"]) is not int or not 0 <= anchor["gene_id"] < len(gene_names)
                or gene_names[anchor["gene_id"]] != anchor["gene"] or anchor["gene"] not in c["members"]
                or anchor["gene"] == name or not math.isfinite(anchor["score"]) or anchor["score"] <= 0
            ):
                raise ValueError("Invalid nominated original sequence anchor")
    if (len(by_family) > 128 or len(selection["selected_profiles"]) != len(by_family)
            or set(selection["selected_profiles"]) != set(by_family)):
        raise ValueError("Complete profile budget/list differs")
    audit_selected = [r["query_id"] for r in selection["audit"] if r["status"] == "selected"]
    if len(audit_selected) != len(by_query) or set(audit_selected) != set(by_query):
        raise ValueError("Selected query audit differs")
    if set(map(int, evidence["measurements"])) != set(by_family):
        raise ValueError("Measurements do not cover the complete competitor set")
    measured = {}
    for family, queries in by_family.items():
        pin = evidence["measurements"][str(family)]
        path = Path(pin["path"]).resolve()
        expected = root / ("builtin_measurements.json" if mode == "builtin_reselect" else f"profile{family:08d}/measurements.json")
        if path != expected:
            raise ValueError("Measurement path is outside its exact arm/profile namespace")
        check(pin)
        checked.append(pin)
        payload = json.loads(path.read_text())
        data = payload[str(family)] if mode == "builtin_reselect" else payload
        expected_status = "builtin_profile_scored" if mode == "builtin_reselect" else "hmmer_profile_scored"
        if (data["status"] != expected_status or data["members"] != groups[family]
                or data["queries"] != sorted(queries) or data["parameters"]["database_genes"] != len(gene_names)
                or set(data["measurements"]) != set(groups[family]) | set(queries)):
            raise ValueError("Backend measurement targets/source/database scope differs")
        if mode != "builtin_reselect":
            for child in [*data["files"], *data["tools"].values()]:
                check(child)
                checked.append(child)
            if expected_tools is not None and data["tools"] != expected_tools:
                raise ValueError("Actual HMMER tools differ from frozen binaries")
            directory = path.parent
            build = [data["tools"]["hmmbuild"]["path"], "--amino", "--hand", "--wnone", "--enone",
                "--seed", "42", "--cpu", "0", "-n", "SelectiveProfile", "-O", str(directory / "hmmer_source.sto"),
                str(directory / "profile.hmm"), str(directory / "source.sto")]
            search = [data["tools"]["hmmsearch"]["path"], "--max", "--seed", "42", "--cpu", "0", "--noali",
                "-T", "-1000000000", "--domT", "-1000000000", "-Z", str(len(gene_names)),
                "--tblout", str(directory / "sequence.tbl"), "--domtblout", str(directory / "domain.tbl"),
                str(directory / "profile.hmm"), str(directory / "targets.faa")]
            if (data["commands"]["build"] != build or data["commands"]["search"] != search
                    or data["parameters"]["sequence_weighting"] != "none"
                    or data["parameters"]["effective_sequence_number"] != len(groups[family])):
                raise ValueError("HMMER construction/search controls differ from frozen protocol")
            targets = {sequence.id: len(sequence) for sequence in SeqIO.parse(path.parent / "targets.faa", "fasta")}
            parsed = read_hmmer_tables(path.parent / "sequence.tbl", path.parent / "domain.tbl",
                data["token_to_gene"], targets, data["model_length"])
            if parsed != data["measurements"]:
                raise ValueError("Raw HMMER tables differ from recorded measurements")
        elif data["minimum_member_score"] != thresholds[family]:
            raise ValueError("Native recalibration differs from original threshold")
        measured[family] = data
    decisions = {}
    for row in evidence["decisions"]:
        q = row["query_id"]
        if q not in by_query or q in decisions or row["query_gene"] != gene_names[q]:
            raise ValueError("Decision query identity differs")
        expected = {"query_id": q, "query_gene": gene_names[q], **decide_query(by_query[q], measured)}
        if row != expected:
            raise ValueError("Decision does not reproduce from complete backend measurements")
        decisions[q] = row
    if set(decisions) != set(by_query):
        raise ValueError("Decisions omit selected queries")
    original = {}
    for row in evidence["original_selected_anchors"]:
        q = row["query_id"]
        if q not in by_query or q in original or row["query_gene"] != gene_names[q]:
            raise ValueError("Invalid original selected anchor identity")
        c = next((c for c in by_query[q]["competitors"] if c["profile_id"] == row["profile_id"]), None)
        if (c is None or c["is_own"] or c["anchor"] != {"gene_id": row["anchor_gene_id"],
            "gene": row["anchor_gene"], "score": row["anchor_score"]}
            or not math.isfinite(row["profile_score"]) or row["profile_score"] < c["original_self_threshold"]):
            raise ValueError("Original anchor differs from nominated source evidence")
        original[q] = row["anchor_gene_id"], row["anchor_score"]
    changed = []
    for q, d in decisions.items():
        after = (original.get(q) if d["action"] == "preserve" else
                 (d["anchor"]["gene_id"], d["anchor"]["score"]) if d["action"] == "recruit" else None)
        if after != original.get(q):
            changed.append(gene_names[q])
    complete = sum(d["action"] != "preserve" for d in decisions.values())
    if (evidence["selected_queries"] != len(by_query) or evidence["completely_measured_queries"] != complete
            or evidence["final_profile_edges"] != graph["counts"]["strict_profile_edges"]
            or evidence["counterfactual_changed_query_anchors"] != sorted(changed)
            or evidence["applied_changed_query_anchors"] != ([] if mode == "hmmer_score_only" else sorted(changed))
            or evidence["measurement_only"] != (mode == "hmmer_score_only")):
        raise ValueError("Actual selected/measured/change/edge counts differ")
    return {"selective_mode": mode, "candidate_digest": selection["candidate_digest"],
            "selected_queries": len(by_query), "selected_profiles": len(by_family),
            "completely_measured_queries": complete, "applied_changed_query_anchors": len(evidence["applied_changed_query_anchors"]),
            "applicability": "completely_measured" if complete else "no_complete_measurements",
            "checked": list({p["path"]: p for p in checked}.values())}
