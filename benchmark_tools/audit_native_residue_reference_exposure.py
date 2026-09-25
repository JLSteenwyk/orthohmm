"""Describe QfO reference exposure of the seven reviewed native O deletions."""

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import subprocess

from benchmark_tools.audit_orthomcl_reference_impact import audit_vgnc, parse_native_report
from benchmark_tools.orthomcl_matrix_to_pairwise import load_species
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.qfo_filter_pairs import load_mapping
from benchmark_tools.reviewed_legacy_database import REVIEW_SHA, validate_content
from benchmark_tools.run_simulation_methods import read_frozen


def select_targets(review, mapping, owners):
    validate_content(dict(status="database_sequence_differences_require_review",
                          content=review["content"]), review)
    rows = []
    for item in review["transformations"]:
        gene = item["id"]
        fields = gene.split("|")
        if len(fields) != 3 or not fields[1] or gene not in owners:
            raise ValueError("Invalid transformation gene identity")
        accession = fields[1]
        number = mapping[accession]
        if type(number) is not int or number < 1:
            raise ValueError("Invalid reference protein number")
        rows.append(dict(accession=accession, gene=gene, protein_number=number,
                         species=owners[gene], transformation=item))
    if len({r["protein_number"] for r in rows}) != len(rows):
        raise ValueError("Non-injective transformation mapping")
    return rows


def summarize(rows, text, vgnc, fas):
    targets = {r["protein_number"] for r in rows}
    families, annotations = parse_native_report(text, targets)
    for family in families:
        family["selected_members"] = family.pop("failed_members")
    records = []
    for row in rows:
        features = fas[row["species"]]["feature"]
        feature = features.get(row["accession"], {})
        records.append(dict(row, **annotations[row["protein_number"]],
            fas_annotation_entry_present=row["accession"] in features,
            fas_feature_types_by_tool={k: len(v) for k, v in feature.items() if isinstance(v, dict)}))
    tree_summary = {}
    for name in ("SwissTrees", "TreeFam-A"):
        selected = [f for f in families if f["benchmark"] == name]
        eligible = [f for f in selected if f["mapped_proteins"] > 5]
        tree_summary[name] = dict(reference_cases=len(selected), eligible_cases=len(eligible),
            eligible_cases_with_selected_members=sum(bool(f["selected_members"]) for f in eligible),
            incident_relation_counts=dict(Counter(r["event"] for f in eligible for r in f["incident_relations"])))
    return dict(status="native_residue_reference_exposure_described", records=records,
        tree_summary=tree_summary, tree_cases_with_selected_members=[f for f in families if f["selected_members"]],
        vgnc=vgnc, selected_proteins=len(rows), exact_sequence_parity=False,
        accuracy_admitted=False, publication_ready=False,
        limitations=["Direct reference exposure is not a counterfactual accuracy effect.",
                     "No query failure, hit coverage or final-group membership is inferred.",
                     "Unexposed proteins may still affect clustering indirectly.",
                     "FAS feature presence does not establish membership in a scored pair."])


def run(root, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    review_path = root / "benchmark_tools/results/qfo_recovery_native_residue_deletions_20260925.json"
    review = read_frozen(review_path, REVIEW_SHA)
    repo = root / "qfo_benchmark/benchmark-webservice"
    reference = repo / "reference_data/2020"
    image = root / "qfo_benchmark/scoring/container_cache/qfobenchmark-darwin-2022.1.img"
    gg = root / "benchmarks/results/qfo_corrected_orthomcl_v1/work/all.gg"
    script = Path(__file__).with_name("orthomcl_reference_impact.drw")
    paths = [review_path, gg, image, script, Path(__file__),
             *[Path(__file__).with_name(n) for n in (
                 "audit_orthomcl_reference_impact.py", "orthomcl_matrix_to_pairwise.py",
                 "prepare_ob_candidate_neighborhood.py", "qfo_filter_pairs.py",
                 "reviewed_legacy_database.py", "run_simulation_methods.py")],
             *[reference / n for n in ("mapping.json.gz", "ReconciledTrees_SwissTrees.drw",
                 "ReconciledTrees_TreeFam-A.drw", "enzymes.drw.gz", "ServerIndexed.db", "vgnc-orthologs.txt.gz")],
             repo / "lib/darwinit"]
    checked = [record(p) for p in paths] + review["checked_records"]
    for item in checked:
        check(item)
    rows = select_targets(review, load_mapping(reference / "mapping.json.gz"), load_species(gg))
    fas_paths = {r["species"]: reference / "fas_annotations" / (r["species"] + ".json") for r in rows}
    checked.extend(record(p) for p in fas_paths.values())
    fas = {s: json.loads(p.read_text()) for s, p in fas_paths.items()}
    output.mkdir(parents=True)
    environment = os.environ.copy()
    targets = "{" + ",".join(str(r["protein_number"]) for r in sorted(rows, key=lambda r: r["protein_number"])) + "}"
    overrides = dict(SINGULARITYENV_QFO_REFSET_PATH=str(reference),
        SINGULARITYENV_DARWIN_ORTHOLOG_BENCHMARK_REPO_PATH=str(repo), SINGULARITYENV_QFO_AUDIT_TARGETS=targets)
    environment.update(overrides)
    command = ["singularity", "exec", str(image), "darwin", "-E"]
    log = output / "native.log"
    with script.open() as source, log.open("x") as stream:
        subprocess.run(command, stdin=source, stdout=stream, stderr=subprocess.STDOUT, env=environment, check=True)
    result = summarize(rows, log.read_text(), audit_vgnc(reference / "vgnc-orthologs.txt.gz",
                       {r["protein_number"] for r in rows}), fas)
    for item in checked:
        check(item)
    result.update(checked_records=checked, native_command=command, environment=overrides, native_log=record(log))
    with (output / "report.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.output.resolve())
