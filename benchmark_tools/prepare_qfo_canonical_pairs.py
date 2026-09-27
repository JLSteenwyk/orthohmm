"""Convert admitted canonical QfO native predictions without changing semantics."""

import argparse
import json
import os
from pathlib import Path
import subprocess

from Bio import SeqIO

from benchmark_tools.readback_qfo_canonical_phylogeny import admit
from benchmark_tools.run_qfo_order_replay import record, save
from benchmark_tools.prepare_qfo_factorial_pairs import write_native_pairs
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.qfo_filter_pairs import load_mapping, filter_pairs
from benchmark_tools.verify_ygob_validation import require_completed_job

PLAN_SHA = "465a0509d7d00640d32cde8411b145ad454cf52627be8469083e37f695fb195d"
SUBMISSION_SHA = "03fed3426661d4c7aa99768c6e67c54f84fea7bea453784b360d66ebd4452637"
AUDIT_PLAN_SHA = "a268529be36550cc0630714a311cceac20b2634ade79d14d4277b46f566a0de9"
RESULT_SHA = "e1409d69b48565ff482a70be83b0e7859ec03be896ebbc015489bd699f6c0ab8"
EXPECTED_PAIRS = 5959535


def verify_execution(execution, plan_record, result_record):
    if execution != dict(status="readback_complete", job_id="22334",
                         plan=plan_record, result=result_record):
        raise ValueError("Wrong canonical scientific audit completion")


def require_counts(count, total, retained):
    if any(type(x) is not int for x in (count, total, retained)) or not count == total == retained == EXPECTED_PAIRS:
        raise ValueError("Canonical pair count mismatch or reference mapping loss")


def check_records(records):
    for item in records:
        if record(item["path"]) != item:
            raise ValueError("Changed conversion dependency: " + item["path"])


def audit_inputs(repo):
    directory = repo / "benchmarks/work/qfo_canonical_phylogeny_20260927"
    accounting = subprocess.check_output(["sacct", "-j", "22334", "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, 22334)
    if scheduler["AllocCPUS"] != "2":
        raise ValueError("Wrong scientific audit allocation")
    audit_path = directory / "readback_plan.json"
    result_path = directory / "readback/result.json"
    audit = read_frozen(audit_path, AUDIT_PLAN_SHA)
    result = read_frozen(result_path, RESULT_SHA)
    execution_path = directory / "readback_execution.json"
    verify_execution(json.loads(execution_path.read_text()), record(audit_path), record(result_path))
    checked = [*audit["checked_records"], *result["reports"], result["admission"],
               *[record(p) for p in (audit_path, result_path, execution_path)]]
    check_records(checked)
    native = admit(repo, directory, 22333, PLAN_SHA, SUBMISSION_SHA)
    if native != json.loads(Path(result["admission"]["path"]).read_text()):
        raise ValueError("Native admission differs from successful scientific audit")
    return directory, native, checked, scheduler


def prepare(repo):
    if not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "2":
        raise ValueError("Require scheduled two-CPU conversion")
    output = repo / "benchmarks/results/qfo_canonical_pairs_20260927"
    if output.exists():
        raise FileExistsError(output)
    directory, native, checked, scheduler = audit_inputs(repo)
    plan = json.loads((directory / "plan.json").read_text())
    phylo = directory / "native/inference/orthohmm_phylogeny"
    manifest_path = phylo / "provenance_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    taxa = {p["filename"]: p["taxon"] for p in manifest["input_proteomes"]}
    if len(taxa) != len(manifest["input_proteomes"]) or set(taxa) != {Path(p["path"]).name for p in plan["fastas"]}:
        raise ValueError("Wrong input proteome inventory")
    owners = {}
    for item in plan["fastas"]:
        for sequence in SeqIO.parse(item["path"], "fasta"):
            if sequence.id in owners:
                raise ValueError("Duplicate input gene")
            owners[sequence.id] = taxa[Path(item["path"]).name]
    if len(owners) != 984137 or len(set(owners.values())) != 78:
        raise ValueError("Wrong complete gene universe")
    environment_path = repo / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(environment_path, ENV_SHA)
    mappings = [r for r in environment["reference_files"] if Path(r["path"]).name == "mapping.json.gz"]
    if len(mappings) != 1:
        raise ValueError("Require one frozen reference mapping")
    pairs_input = record(phylo / "orthohmm_pairwise_orthologs.tsv")
    if pairs_input not in native["outputs"]:
        raise ValueError("Pair input not admitted")
    checked.extend([pairs_input, record(manifest_path), record(environment_path), mappings[0], record(__file__)])
    check_records(checked)
    output.mkdir(parents=True)
    report = dict(status="preparing", participant="ohmm_qfo_canonical_20260927",
        job_id=os.environ["SLURM_JOB_ID"], audit_scheduler=scheduler,
        semantics="native phylogenetically inferred pairs", native_input=pairs_input,
        mapping=mappings[0], input_fastas=plan["fastas"], checked_records=checked,
        accuracy_evaluated=False, publication_ready=False)
    save(output / "preflight.json", report)
    try:
        pairs, filtered = output / "pairs.partial.tsv", output / "pairs.qfo.partial.tsv"
        count = write_native_pairs(Path(pairs_input["path"]), pairs, owners, EXPECTED_PAIRS)
        total, retained = filter_pairs(pairs, filtered, load_mapping(Path(mappings[0]["path"])))
        require_counts(count, total, retained)
        check_records(checked)
        final, mapped = output / "pairs.tsv", output / "pairs.qfo.tsv"
        pairs.rename(final)
        filtered.rename(mapped)
        report.update(status="canonical_native_pairs_prepared_unscored", pairs=record(final),
            filtered_pairs=record(mapped), total_pairs=total, retained_pairs=retained, removed_mapping_pairs=0)
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save(output / "results.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.repo.resolve())
    print(json.dumps(dict(status=result["status"], pairs=result["filtered_pairs"])))
