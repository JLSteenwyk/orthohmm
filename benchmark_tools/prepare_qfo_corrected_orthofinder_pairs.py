"""Convert admitted corrected OrthoFinder native pairs or MCL diagnostic groups."""

import argparse
from collections import Counter, defaultdict
from itertools import combinations, product
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_high_sensitivity import PLAN_SHA
from benchmark_tools.audit_orthofinder_pair_tables import read_table
from benchmark_tools.orthofinder_to_pairwise import _accession
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.qfo_filter_pairs import load_mapping, filter_pairs
from benchmark_tools.report_ygob_validation import read_checkpoint
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.score_ygob_groups import membership
from benchmark_tools.validate_scaling_outputs import input_universe
from benchmark_tools.verify_ygob_validation import require_completed_job

ADMITTER = "33e2310d1095f64e77cf068ef2a6fccd4b32a8d4"
METHODS = ("orthofinder_full", "orthofinder_sequence_only")
SEMANTICS = {METHODS[0]: "native phylogenetically inferred pairs",
             METHODS[1]: "cross-species pre-phylogenetic MCL group-derived clique pairs (diagnostic)"}


def validate_admission(admission):
    if (admission["status"] != "corrected_orthofinder_native_evidence_admitted"
            or admission["accuracy_evaluated"] is not False or admission["publication_ready"] is not False):
        raise ValueError("Require corrected native OrthoFinder admission")
    content, scheduler = admission["content"], admission["scheduler"]
    if (content["genes"] != 984137 or content["species"] != 78
            or scheduler["State"] != "COMPLETED" or scheduler["ExitCode"] != "0:0"
            or scheduler["AllocCPUS"] != "32" or scheduler["NodeList"] != "bizon"):
        raise ValueError("Wrong corrected input universe or inference allocation")
    pairs = content["native_pairs"]
    if any(type(pairs[k]) is not int or pairs[k] < 0 for k in (
            "distinct_pairs", "converter_emitted_pairs", "converter_duplicate_pairs")):
        raise ValueError("Invalid native pair counts")
    if (pairs["status"] != "native_orthofinder_pair_tables_verified"
            or pairs["mcl_membership_checked"] is not True or pairs["directed_tables"] != 6006
            or pairs["species"] != 78 or pairs["input_genes"] != 984137
            or pairs["distinct_pairs"] <= 0
            or pairs["converter_emitted_pairs"] != pairs["distinct_pairs"] + pairs["converter_duplicate_pairs"]):
        raise ValueError("Invalid native pair audit")
    return content


def accessions(owners):
    result = {g: _accession(g) for g in owners}
    if any(not a for a in result.values()) or len(set(result.values())) != len(result):
        raise ValueError("Accession conversion is not injective")
    return result


def write_native(content, owners, groups, stream):
    aliases, member = accessions(owners), membership(groups)
    species = sorted(set(owners.values()))
    rows = content["native_pairs"]["species_pairs"]
    if [r["species"] for r in rows] != [list(p) for p in combinations(species, 2)]:
        raise ValueError("Changed native species-pair inventory/order")
    total, duplicates = 0, 0
    for row in rows:
        a, b = row["species"]
        path = Path(content["results_directory"]) / "Orthologues" / f"Orthologues_{a}" / f"{a}__v__{b}.tsv"
        pairs, stats = read_table(path, a, b, owners, member)
        if stats != row["forward"]:
            raise ValueError("Native relation counts changed since admission")
        for pair in sorted(pairs):
            x, y = sorted(aliases[g] for g in pair)
            stream.write(f"{x}\t{y}\n")
            total += 1
        duplicates += stats["duplicate_relations"]
    if (total != content["native_pairs"]["distinct_pairs"]
            or duplicates != content["native_pairs"]["converter_duplicate_pairs"]):
        raise ValueError("Native distinct/duplicate pair counts differ")
    return total, duplicates


def write_groups(groups, owners, stream):
    aliases = accessions(owners)
    if set(membership(groups)) != set(owners):
        raise ValueError("MCL partition does not cover the exact input universe")
    total, expected = 0, 0
    for genes in groups.values():
        counts = Counter(owners[g] for g in genes)
        expected += (len(genes) ** 2 - sum(n * n for n in counts.values())) // 2
        buckets = defaultdict(list)
        for gene in genes:
            buckets[owners[gene]].append(aliases[gene])
        for a, b in combinations(sorted(buckets), 2):
            for pair in product(buckets[a], buckets[b]):
                x, y = sorted(pair)
                stream.write(f"{x}\t{y}\n")
                total += 1
    if total != expected:
        raise ValueError("MCL pair expansion differs from independent species-size formula")
    return total, 0


def prepare(root, method, admission_path, admission_sha, admission_job):
    if method not in METHODS:
        raise ValueError("Unknown OrthoFinder output semantics")
    if os.environ.get("SLURM_CPUS_PER_TASK") != "2" or not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("Require scheduled two-CPU conversion")
    accounting = subprocess.check_output(["sacct", "-j", str(admission_job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, admission_job)
    if scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "2":
        raise ValueError("Wrong native admission allocation")
    admission = read_frozen(admission_path, admission_sha)
    content = validate_admission(admission)
    executor = root / "benchmarks/work/publication_qfo_corrected_orthofinder_admission_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != ADMITTER:
        raise ValueError("Native admission executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    if admission["source"] != record(executor / "benchmark_tools/admit_qfo_corrected_orthofinder.py"):
        raise ValueError("Wrong native admission source")
    plan_path = root / "benchmark_tools/results/qfo_corrected_primary_commands_20260918.json"
    plan = read_frozen(plan_path, PLAN_SHA)
    if record(plan_path) not in admission["checked_records"]:
        raise ValueError("Native admission not bound to corrected primary plan")
    inputs = [r for r in plan["inputs"] if Path(r["path"]).parent == Path(plan["input_directory"])
              and Path(r["path"]).suffix == ".fasta"]
    env_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(env_path, ENV_SHA)
    mappings = [r for r in environment["reference_files"] if Path(r["path"]).name == "mapping.json.gz"]
    if len(mappings) != 1:
        raise ValueError("Require one frozen QfO mapping")
    mapping = mappings[0]
    checked = [record(admission_path), admission["source"], *admission["checked_records"],
               record(plan_path), record(env_path), mapping, record(__file__),
               *[record(Path(__file__).with_name(n)) for n in (
                   "audit_orthofinder_pair_tables.py", "orthofinder_to_pairwise.py",
                   "report_ygob_validation.py", "orthofinder_mcl_to_orthogroups.py",
                   "score_ygob_groups.py", "validate_scaling_outputs.py", "qfo_filter_pairs.py")]]
    for item in checked:
        check(item)
    owners, _ = input_universe({"inputs": inputs, "proteins": 984137, "proteomes": 78})
    groups = read_checkpoint(Path(content["checkpoint"]["path"]), Path(content["sequence_ids"]["path"]), owners)
    if len(groups) != content["checkpoint_groups"]:
        raise ValueError("Checkpoint group count changed")
    directory = root / "benchmarks/results/qfo_corrected_comparator_pairs_v1" / method
    directory.mkdir(parents=True, exist_ok=False)
    report = {"status": "preparing", "source": record(__file__), "method": method,
              "participant": "qfo_corrected_" + method, "admission": record(admission_path),
              "admission_scheduler": scheduler, "mapping": mapping, "checked_records": checked,
              "accuracy_evaluated": False, "publication_ready": False,
              "job_id": os.environ["SLURM_JOB_ID"], "started_epoch": time.time(), "semantics": SEMANTICS[method]}
    try:
        partial, filtered_partial = directory / "pairs.partial.tsv", directory / "pairs.qfo.partial.tsv"
        with partial.open("x") as stream:
            total, duplicates = (write_native(content, owners, groups, stream) if method == METHODS[0]
                                 else write_groups(groups, owners, stream))
        observed, retained = filter_pairs(partial, filtered_partial, load_mapping(Path(mapping["path"])))
        if not 0 < total == observed == retained:
            raise ValueError("Pair count mismatch or unexpected corrected-reference mapping loss")
        for item in checked:
            check(item)
        pairs, filtered = directory / "pairs.tsv", directory / "pairs.qfo.tsv"
        partial.rename(pairs)
        filtered_partial.rename(filtered)
        report.update(status="corrected_orthofinder_pairs_prepared_unscored", pairs=record(pairs),
                      filtered_pairs=record(filtered), total_pairs=total, retained_pairs=retained,
                      removed_mapping_pairs=0, native_duplicate_relations=duplicates)
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
    for flag in ("root", "admission"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--method", choices=METHODS, required=True)
    parser.add_argument("--admission-sha256", required=True)
    parser.add_argument("--admission-job", type=int, required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.method, args.admission.resolve(), args.admission_sha256, args.admission_job)
