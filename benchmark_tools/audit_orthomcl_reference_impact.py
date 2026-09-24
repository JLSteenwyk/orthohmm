"""Measure direct reference exposure of audited OrthoMCL failed queries."""

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.qfo_filter_pairs import load_mapping


def failed_query_targets(audit, mapping, allow_grouped=False):
    """Keep failed search queries distinct from absence in the final partition."""
    failed = []
    genes, accessions = set(), set()
    for item in audit["records"]:
        if type(item["query_failed"]) is not bool:
            raise ValueError("Invalid query failure flag")
        if not item["query_failed"]:
            continue
        gene, accession = item["gene"], item["accession"]
        if not gene or not accession or gene in genes or accession in accessions:
            raise ValueError("Missing/duplicate failed-query identity")
        genes.add(gene)
        accessions.add(accession)
        line, size = item["final_group_line"], item["final_group_size"]
        if (type(size) is not int or (line is None and size != 0)
                or (line is not None and (type(line) is not int or line < 1 or size < 1))):
            raise ValueError("Inconsistent failed-query final-group membership")
        if line is not None and not allow_grouped:
            raise ValueError("Audit does not establish missing failed-query predictions; opt in to grouped exposure")
        failed.append(item)
    if type(audit["failed_queries"]) is not int or len(failed) != audit["failed_queries"]:
        raise ValueError("Failed-query count differs from records")
    targets = [mapping[item["accession"]] for item in failed]
    if (not targets or any(type(p) is not int or p < 1 for p in targets)
            or len(set(targets)) != len(targets)):
        raise ValueError("Empty, invalid or non-injective mapped failed-query targets")
    return failed, set(targets)


def parse_native_report(text, targets):
    if re.search(r"^\s*(?:Error|ERROR|error)\b", text, re.MULTILINE):
        raise ValueError("Native interpreter reported an error")
    families = {}
    annotations = {}
    complete = False
    for line in text.splitlines():
        if not line.startswith("QFOAUDIT\t"):
            continue
        fields = line.split("\t")[1:]
        kind = fields[0]
        if kind == "FAMILY" and len(fields) == 4:
            key = (fields[1], fields[2])
            if key in families:
                raise ValueError("Duplicate native family")
            families[key] = {"benchmark": key[0], "family": key[1],
                             "mapped_proteins": int(fields[3]), "failed_members": [],
                             "incident_relations": []}
        elif kind == "MEMBER" and len(fields) == 4:
            gene = int(fields[3])
            if gene not in targets:
                raise ValueError("Unexpected native member")
            families[(fields[1], fields[2])]["failed_members"].append(gene)
        elif kind == "RELATION" and len(fields) == 6:
            pair = [int(fields[3]), int(fields[4])]
            if pair[0] >= pair[1] or not set(pair) & targets or fields[5] not in {"S", "D", "OTHER"}:
                raise ValueError("Invalid native relation")
            families[(fields[1], fields[2])]["incident_relations"].append({"pair": pair, "event": fields[5]})
        elif kind == "ANNOTATION" and len(fields) == 4:
            gene = int(fields[1])
            if gene in annotations or gene not in targets:
                raise ValueError("Invalid annotation ID")
            values = [int(v) for v in fields[2:]]
            if min(values) < 0:
                raise ValueError("Invalid annotation counts")
            annotations[gene] = dict(zip(("ec_terms", "experimental_go_terms"), values))
        elif kind == "COMPLETE" and len(fields) == 2:
            if complete or int(fields[1]) != len(targets):
                raise ValueError("Invalid completion marker")
            complete = True
        else:
            raise ValueError("Unrecognized native report record")
    if not complete or set(annotations) != targets:
        raise ValueError("Incomplete native report")
    if {key[0] for key in families} != {"SwissTrees", "TreeFam-A"}:
        raise ValueError("Missing reference-tree benchmark")
    return list(families.values()), annotations


def audit_vgnc(path, targets):
    pairs = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            left, right, family = line.rstrip("\n").split("\t")
            pair = tuple(sorted((int(left), int(right))))
            if pair[0] == pair[1] or not family:
                raise ValueError("Invalid VGNC reference pair")
            if pair in pairs and pairs[pair] != family:
                raise ValueError("Conflicting VGNC pair")
            pairs[pair] = family
    incident = [{"pair": list(pair), "family": family}
                for pair, family in sorted(pairs.items()) if set(pair) & targets]
    return {"reference_pairs": len(pairs), "incident_pairs": incident,
            "incident_pair_fraction": len(incident) / len(pairs) if pairs else 0.0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--benchmark-repo", type=Path, required=True)
    parser.add_argument("--darwin-image", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-grouped-failed-queries", action="store_true",
                        help="Measure reference exposure even when failed queries occur in final groups")
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Refusing to overwrite existing impact audit")
    audit = json.loads(args.audit.read_text())
    mapping = load_mapping(args.reference / "mapping.json.gz")
    failed, targets = failed_query_targets(audit, mapping, args.allow_grouped_failed_queries)
    args.output.mkdir(parents=True)
    script = Path(__file__).with_name("orthomcl_reference_impact.drw")
    environment = os.environ.copy()
    environment.update(
        SINGULARITYENV_QFO_REFSET_PATH=str(args.reference.resolve()),
        SINGULARITYENV_DARWIN_ORTHOLOG_BENCHMARK_REPO_PATH=str(args.benchmark_repo.resolve()),
        SINGULARITYENV_QFO_AUDIT_TARGETS="{" + ",".join(map(str, sorted(targets))) + "}",
    )
    command = ["singularity", "exec", str(args.darwin_image.resolve()), "darwin", "-E"]
    with script.open() as source, (args.output / "native.log").open("w") as log:
        subprocess.run(command, stdin=source, stdout=log, stderr=subprocess.STDOUT,
                       env=environment, check=True)
    families, annotations = parse_native_report((args.output / "native.log").read_text(), targets)
    vgnc = audit_vgnc(args.reference / "vgnc-orthologs.txt.gz", targets)
    fas_presence = {}
    fas_features = {}
    fas_files = []
    for species in sorted({r["species"] for r in failed}):
        path = args.reference / "fas_annotations" / f"{species}.json"
        if not path.exists():
            raise ValueError(f"Missing FAS annotation input: {path}")
        data = json.loads(path.read_text())
        fas_presence.update({r["accession"]: r["accession"] in data["feature"]
                             for r in failed if r["species"] == species})
        for record in failed:
            if record["species"] == species:
                feature = data["feature"].get(record["accession"], {})
                fas_features[record["accession"]] = {
                    tool: len(features) for tool, features in feature.items()
                    if isinstance(features, dict)
                }
        fas_files.append(file_provenance(path))
    records = [{"accession": r["accession"], "protein_number": mapping[r["accession"]],
                "length": r["length"], "species": r["species"],
                "final_group_line": r["final_group_line"], "final_group_size": r["final_group_size"],
                **annotations[mapping[r["accession"]]],
                "fas_annotation_entry_present": fas_presence[r["accession"]],
                "fas_feature_types_by_tool": fas_features[r["accession"]]}
               for r in failed]
    tree_summary = {}
    for name in ("SwissTrees", "TreeFam-A"):
        selected = [f for f in families if f["benchmark"] == name]
        eligible = [f for f in selected if f["mapped_proteins"] > 5]
        tree_summary[name] = {
            "reference_cases": len(selected), "eligible_cases": len(eligible),
            "eligible_cases_with_failed_members": sum(bool(f["failed_members"]) for f in eligible),
            "incident_relation_counts": dict(Counter(r["event"] for f in eligible for r in f["incident_relations"])),
        }
    reference_paths = [args.reference / name for name in (
        "mapping.json.gz", "ReconciledTrees_SwissTrees.drw", "ReconciledTrees_TreeFam-A.drw",
        "enzymes.drw.gz", "ServerIndexed.db", "vgnc-orthologs.txt.gz")]
    result = {
        "schema_version": 1, "generated_at": datetime.now(timezone.utc).isoformat(),
        "command": [sys.executable, *sys.argv], "native_command": command,
        "failed_queries": len(failed), "unique_mapped_targets": len(targets),
        "failed_query_group_coverage": {
            "grouped": sum(r["final_group_line"] is not None for r in failed),
            "ungrouped": sum(r["final_group_line"] is None for r in failed),
            "grouped_exposure_enabled": args.allow_grouped_failed_queries},
        "tree_summary": tree_summary, "vgnc": vgnc, "records": records,
        "tree_cases_with_failed_members": [f for f in families if f["failed_members"]],
        "annotation_summary": {
            "failed_queries_with_ec": sum(r["ec_terms"] > 0 for r in records),
            "failed_queries_with_experimental_go": sum(r["experimental_go_terms"] > 0 for r in records),
            "failed_queries_with_fas_entry": sum(r["fas_annotation_entry_present"] for r in records),
            "failed_queries_with_any_fas_feature": sum(any(r["fas_feature_types_by_tool"].values()) for r in records)},
        "inputs": [file_provenance(p) for p in [args.audit, script, Path(__file__), args.darwin_image,
                    args.benchmark_repo / "lib/darwinit", args.benchmark_repo / "RefPhyloTest.drw",
                    args.benchmark_repo / "GoTest.drw", args.benchmark_repo / "EcTest.drw",
                    args.benchmark_repo / "vgnc_benchmark.py", args.benchmark_repo / "fas_benchmark.py",
                    args.benchmark_repo / "nextflow.config",
                    *reference_paths]],
        "fas_inputs": fas_files, "native_log": file_provenance(args.output / "native.log"),
        "limitations": [
            "This is direct reference exposure, not a counterfactual rerun or bound on clustering changes.",
            "A failed outgoing search can coexist with incoming hits or final-group membership; neither proves repaired orthology.",
            "Incident reference pairs need not become correct predictions after repairing search.",
            "GO counts use the six experimental evidence codes selected by the benchmark configuration.",
            "FAS entry presence is not proof of a scored pair or a nonempty feature architecture.",
            "FAS precomputed-pair availability is not assessed here.",
            "Native tree cases are not necessarily individual biological families; TreeFam-A is a pooled case.",
        ],
    }
    (args.output / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: result[k] for k in ("failed_queries", "tree_summary", "annotation_summary", "vgnc")}, indent=2))


if __name__ == "__main__":
    main()
