"""Expand admitted native OrthoMCL final groups into QfO cross-species cliques."""

import argparse
from collections import defaultdict
from itertools import combinations, product
import json
import os
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_orthomcl import completed, frozen, unique_records
from benchmark_tools.audit_orthomcl_native_groups import audit
from benchmark_tools.normalize_three_kingdoms_orthogroups import iter_orthomcl
from benchmark_tools.orthomcl_matrix_to_pairwise import load_species, _accession
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.qfo_filter_pairs import load_mapping, filter_pairs
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_orthomcl_python_runtime import verify_runtime

ADMITTER = "fe08c7d99f266f2705e32c85294adfec676c008f"
SEMANTICS = "cross_species_final_group_cliques"


def validate_admission(admission):
    if (admission["status"] != "corrected_orthomcl_native_outputs_admitted"
            or admission["accuracy_admitted"] is not False or admission["publication_ready"] is not False
            or admission["pair_semantics"] != SEMANTICS):
        raise ValueError("Require native final-group admission, not graph edges")
    content = admission["content"]
    required = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "180", "ReqMem": "900G"}
    if (content["input_proteins"] != 984137 or content["input_species"] != 78
            or any(admission["scheduler"].get(k) != v for k, v in required.items())):
        raise ValueError("Wrong corrected input scope or native scheduler")
    for key in ("cross_species_clique_pairs", "final_groups", "grouped_proteins", "ungrouped_input_proteins"):
        if type(content[key]) is not int or content[key] < 0:
            raise ValueError("Invalid admitted pair/group counts")
    if content["grouped_proteins"] + content["ungrouped_input_proteins"] != 984137:
        raise ValueError("Inconsistent admitted group coverage")
    return content


def write_pairs(groups, owners, stream):
    """Stream disjoint group cliques; global accession collisions are errors."""
    aliases = {gene: _accession(gene) for gene in owners}
    if any(not value for value in aliases.values()) or len(set(aliases.values())) != len(aliases):
        raise ValueError("Non-injective accession mapping")
    seen, labels = set(), set()
    total = expected = 0
    for label, genes in iter_orthomcl(groups):
        members = set(genes)
        if (label in labels or len(members) != len(genes) or len(genes) < 2
                or not members <= owners.keys() or members & seen):
            raise ValueError("Unknown, singleton or duplicate final-group membership")
        labels.add(label)
        seen.update(members)
        buckets = defaultdict(list)
        for gene in genes:
            buckets[owners[gene]].append(aliases[gene])
        expected += (len(genes)**2 - sum(len(values)**2 for values in buckets.values())) // 2
        for a, b in combinations(sorted(buckets), 2):
            for pair in product(sorted(buckets[a]), sorted(buckets[b])):
                x, y = sorted(pair)
                stream.write(f"{x}\t{y}\n")
                total += 1
    if total != expected:
        raise ValueError("Clique expansion differs from species-size formula")
    return {"total_pairs": total, "final_groups": len(labels), "grouped_proteins": len(seen),
            "ungrouped_input_proteins": len(owners) - len(seen)}


def prepare(root, admission_path, admission_sha, admission_job):
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "2"
            or os.environ.get("SLURM_MEM_PER_NODE") != "65536" or os.uname().nodename != "bizon"):
        raise ValueError("Require scheduled 2-CPU/64-GiB conversion on bizon")
    runtime = verify_runtime(root)
    scheduler, accounting = completed(admission_job, 2, "64G")
    expected_path = root / "benchmarks/work/qfo_corrected_orthomcl_admission_20260918/report.json"
    if admission_path != expected_path:
        raise ValueError("Wrong native admission path")
    admission = read_frozen(admission_path, admission_sha)
    content = validate_admission(admission)
    executor = root / "benchmarks/work/publication_qfo_corrected_orthomcl_admission_v1"
    frozen(executor, ADMITTER)
    if admission["source"] != record(executor / "benchmark_tools/admit_qfo_corrected_orthomcl.py"):
        raise ValueError("Wrong native admission source")
    env_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(env_path, ENV_SHA)
    mappings = [r for r in environment["reference_files"] if Path(r["path"]).name == "mapping.json.gz"]
    if len(mappings) != 1:
        raise ValueError("Require one frozen reference mapping")
    mapping = mappings[0]
    checked = unique_records([record(admission_path), admission["source"], *admission["checked_records"],
        *admission["outputs"], admission["native_groups"], record(env_path), mapping,
        *[record(Path(__file__).with_name(name)) for name in (
            "prepare_qfo_corrected_orthomcl_pairs.py", "admit_qfo_corrected_orthomcl.py",
            "audit_orthomcl_native_groups.py", "normalize_three_kingdoms_orthogroups.py",
            "orthomcl_matrix_to_pairwise.py", "qfo_filter_pairs.py")]])
    for item in checked:
        check(item)
    base = root / "benchmarks/results/qfo_corrected_orthomcl_v1"
    native = Path(admission["native_groups"]["path"]).parent
    if native.parent != base / "native_tool" or admission["native_groups"]["path"] != str(native / "all_orthomcl.out"):
        raise ValueError("Unexpected native final-group path")
    gg = base / "inference_execution/inputs/all.gg"
    owners = load_species(gg)
    if len(owners) != 984137 or len(set(owners.values())) != 78:
        raise ValueError("Changed corrected GG scope")
    directory = root / "benchmarks/results/qfo_corrected_comparator_pairs_v1/orthomcl"
    directory.mkdir(parents=True, exist_ok=False)
    report = {"status": "preparing", "source": record(__file__), "method": "orthomcl",
              "participant": "qfo_corrected_orthomcl", "admission": record(admission_path),
              "admission_scheduler": scheduler, "admission_accounting": accounting,
              "mapping": mapping, "checked_records": checked, "runtime_before": runtime,
              "query_coverage": admission["query_coverage"], "accuracy_evaluated": False,
              "publication_ready": False, "job_id": os.environ["SLURM_JOB_ID"],
              "started_epoch": time.time(), "semantics": SEMANTICS}
    try:
        validation = audit(native / "all_orthomcl.out", native / "tmp/all_ortho.mcl",
                           native / "tmp/all_ortho.idx", gg, directory / "groups.json")
        if validation["content"] != content:
            raise ValueError("Native final-group audit differs from admission")
        partial, mapped_partial = directory / "pairs.partial.tsv", directory / "pairs.qfo.partial.tsv"
        with partial.open("x") as stream:
            counts = write_pairs(native / "all_orthomcl.out", owners, stream)
        expected = {key: content[key] for key in ("final_groups", "grouped_proteins", "ungrouped_input_proteins")}
        expected["total_pairs"] = content["cross_species_clique_pairs"]
        if counts != expected:
            raise ValueError("Conversion counts differ from admission")
        observed, retained = filter_pairs(partial, mapped_partial, load_mapping(Path(mapping["path"])))
        if observed != retained or observed != counts["total_pairs"]:
            raise ValueError("Unexpected mapping loss or pair count mismatch")
        for item in checked:
            check(item)
        for item in validation["checked_records"]:
            check(item)
        runtime_after = verify_runtime(root)
        pairs, filtered = directory / "pairs.tsv", directory / "pairs.qfo.tsv"
        partial.rename(pairs)
        mapped_partial.rename(filtered)
        report.update(status="corrected_orthomcl_pairs_prepared_unscored", content=counts,
                      pairs=record(pairs), filtered_pairs=record(filtered), group_audit=record(directory / "groups.json"),
                      total_pairs=observed, retained_pairs=retained, removed_mapping_pairs=0,
                      native_duplicate_relations=0, runtime_after=runtime_after)
        report["limitations"] = ["Cross-species final-group cliques, not native phylogenetic pairs or pre-clustering edges.",
                                  "Ungrouped proteins are not silently added or counted as paired predictions.",
                                  "Source query failures remain retained; no repaired hits are implied.",
                                  "Conversion is unscored and requires terminal verification before QfO assessment."]
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["finished_epoch"] = time.time()
        with (directory / "results.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "admission"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    parser.add_argument("--admission-job", type=int, required=True)
    args = parser.parse_args()
    print(json.dumps({"status": prepare(args.root.resolve(), args.admission.resolve(),
          args.admission_sha256, args.admission_job)["status"]}))
