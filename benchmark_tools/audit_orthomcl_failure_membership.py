"""Join admitted recovered search failures to native final-group membership."""

import argparse
from collections import Counter
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_orthomcl_native_groups import HEADER, partition, raw_index, validate
from benchmark_tools.orthomcl_matrix_to_pairwise import load_species
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.prepare_blast_recovery_bpo import verify_admission


def describe(search, groups, mcl, index, gg):
    content = validate(groups, mcl, index, gg)
    owners = load_species(gg)
    genes = raw_index(index, owners)
    assignments, sizes = partition(mcl, len(genes))
    positions = {gene: i for i, gene in enumerate(genes)}
    lines = {}
    with groups.open() as stream:
        for number, line in enumerate(stream, 1):
            match = HEADER.fullmatch(line.rstrip("\r\n"))
            lines[int(match[1])] = number
    if search["input_proteins"] != len(owners):
        raise ValueError("Search and native input universes differ")
    rows, seen = [], set()
    for diagnostic in search["diagnostics"]:
        gene = diagnostic["gene"]
        if gene not in owners or gene in seen:
            raise ValueError("Unknown or repeated diagnostic gene")
        seen.add(gene)
        flags = ("query_failed", "has_query_hits", "has_subject_hits", "has_self_hit")
        if (any(type(diagnostic[key]) is not bool for key in flags)
                or type(diagnostic["length"]) is not int or diagnostic["length"] <= 0):
            raise ValueError("Invalid diagnostic flags or length")
        failed = any(m["category"].endswith("failure") for m in diagnostic["messages"])
        if failed != diagnostic["query_failed"] or (failed and diagnostic["has_query_hits"]):
            raise ValueError("Inconsistent failed-query evidence")
        if diagnostic["has_self_hit"] and not (diagnostic["has_query_hits"] and diagnostic["has_subject_hits"]):
            raise ValueError("Self hit contradicts directional hit flags")
        cid = assignments[positions[gene]] if gene in positions else None
        size = sizes[cid] if cid is not None else 0
        state = "absent_from_index" if cid is None else "mcl_singleton" if size == 1 else "final_group"
        fields = gene.split("|")
        accession = fields[1] if len(fields) == 3 else gene
        rows.append(dict(diagnostic, accession=accession, species=owners[gene],
            membership_state=state, native_cluster_id=cid,
            final_group_line=lines[cid] if size > 1 else None,
            final_group_size=size if size > 1 else 0))
    failed = [r for r in rows if r["query_failed"]]
    expected = dict(failed_queries=len(failed), failed_queries_with_query_hits=0,
                    failed_queries_with_subject_hits=sum(r["has_subject_hits"] for r in failed))
    if any(type(search[k]) is not int or search[k] != v for k, v in expected.items()):
        raise ValueError("Failure totals disagree with diagnostic records")
    return dict(input_proteins=len(owners), input_species=len(set(owners.values())),
        native_content=content, **expected, failed_queries_in_final_groups=sum(r["final_group_line"] is not None for r in failed),
        failed_queries_by_membership=dict(Counter(r["membership_state"] for r in failed)),
        records=sorted(rows, key=lambda r: r["gene"]))


def audit(search_path, search_sha, native_path, native_sha, native_representation_root=None):
    search = read_frozen(search_path, search_sha)
    native = read_frozen(native_path, native_sha)
    expected_status = "recovered_orthomcl_search_evidence_verified"
    representation_records = []
    if native_representation_root is not None:
        expected_status = "recovered_search_native_representation_verified"
        verified, _, representation_records, _, _ = verify_admission(
            native_representation_root, search_path, search_sha)
        if verified != search:
            raise ValueError("Search changed during representation verification")
    if (search["status"] != expected_status or search["search_admitted"] is not True
            or native["status"] != "recovered_orthomcl_native_outputs_admitted"
            or native["conversion_authorized"] is not True
            or native["query_coverage"] != search["query_coverage"]):
        raise ValueError("Require matching recovered search and native admissions")
    group_audits = [r for r in native["outputs"] if Path(r["path"]).name == "groups.json"]
    if len(group_audits) != 1:
        raise ValueError("Missing or ambiguous native group audit")
    group_record = group_audits[0]
    group_audit = read_frozen(Path(group_record["path"]), group_record["sha256"])
    if (group_audit["status"] != "native_final_groups_match_mcl_partition"
            or group_audit["content"] != native["content"]):
        raise ValueError("Group audit disagrees with native admission")
    selected = []
    for name in ("all_orthomcl.out", "all_ortho.mcl", "all_ortho.idx", "all.gg"):
        matches = [r for r in group_audit["checked_records"] if Path(r["path"]).name == name]
        if len(matches) != 1:
            raise ValueError("Missing or ambiguous group input")
        selected.append(matches[0])
    if selected[0] != native["native_groups"]:
        raise ValueError("Final groups differ from admitted output")
    checked = [record(search_path), record(native_path), group_record,
               *group_audit["checked_records"], *representation_records]
    for item in checked:
        check(item)
    result = describe(search["query_coverage"], *[Path(r["path"]) for r in selected])
    if result["native_content"] != native["content"]:
        raise ValueError("Reconstructed partition differs from admission")
    for item in checked:
        check(item)
    result.update(status="recovered_failure_membership_content_verified", checked_records=checked,
        source=record(__file__), helpers=[record(Path(__file__).with_name(n)) for n in (
            "audit_orthomcl_native_groups.py", "orthomcl_matrix_to_pairwise.py",
            "prepare_ob_candidate_neighborhood.py", "run_simulation_methods.py")],
        accuracy_admitted=False, publication_ready=False,
        limitations=["Relies on pinned upstream search/native admissions; does not rerun search or validate scheduler execution.",
            "Incoming hits and final-group membership do not repair failed outgoing search or prove orthology.",
            "Native singleton clusters are not added to final output; unindexed proteins remain distinct.",
            "No reference exposure, mapping validation or counterfactual accuracy effect is measured here."])
    if native_representation_root is not None:
        result["database_representation"] = search["database_representation"]
        result["limitations"].extend(search["database_representation"]["limitations"])
        result["helpers"].extend(record(Path(__file__).with_name(n)) for n in (
            "prepare_blast_recovery_bpo.py", "reviewed_legacy_database.py"))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("search", "native", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("search", "native"):
        parser.add_argument("--" + name + "-sha256", required=True)
    parser.add_argument("--native-representation-root", type=Path,
                        help="Explicitly verify the reviewed legacy-residue search contract")
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    root = args.native_representation_root.resolve() if args.native_representation_root else None
    result = audit(args.search.resolve(), args.search_sha256, args.native.resolve(), args.native_sha256, root)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
