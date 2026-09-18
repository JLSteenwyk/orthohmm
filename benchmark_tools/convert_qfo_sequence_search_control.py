"""Convert an admitted corrected QfO search panel into numeric checkpoints."""

import argparse
import json
import os
from pathlib import Path
import sqlite3
import subprocess
import sys

from Bio import SeqIO
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_sequence_search_control import validate_execution, EXECUTOR
from benchmark_tools.convert_sequence_search_control import initialize, ingest, parse_hits, write_variant
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_qfo_sequence_search_control import verify_plan
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_ygob_validation import require_completed_job

ADMITTER = "6cf2e92074733cab036c1dd1674467bd37ce2313"


def validate_metadata(inputs, metadata, expected_genes=984137, expected_species=78):
    seen, species = set(), set()
    for item in inputs:
        path = Path(item["path"])
        check(item)
        species.add(path.name)
        for sequence in SeqIO.parse(path, "fasta"):
            value = metadata.get(sequence.id)
            if (sequence.id in seen or not len(sequence.seq) or not isinstance(value, dict)
                    or set(value) != {"length", "species"} or type(value["length"]) is not int
                    or value["length"] != len(sequence.seq) or value["species"] != path.name):
                raise ValueError("Metadata differs from exact corrected FASTA IDs, lengths or ownership")
            seen.add(sequence.id)
    if seen != set(metadata) or len(seen) != expected_genes or len(species) != expected_species:
        raise ValueError("Incomplete corrected input universe")


def convert_hits(plan, metadata, output):
    names = sorted(metadata)
    indices = {name: i for i, name in enumerate(names)}
    owners = {name: i for i, name in enumerate(sorted({r["species"] for r in metadata.values()}))}
    species = np.asarray([owners[metadata[name]["species"]] for name in names], dtype=np.int32)
    with sqlite3.connect(output / "hits.sqlite") as database:
        database.execute("PRAGMA temp_store=FILE")
        database.execute("PRAGMA cache_size=-131072")
        initialize(database)
        for target in plan["searches"]:
            owner = Path(target["target_fasta"]["path"]).name
            ingest(database, parse_hits(Path(target["output"]), metadata, indices, owner, owners[owner]))
        if database.execute("PRAGMA quick_check").fetchall() != [("ok",)]:
            raise ValueError("SQLite integrity check failed")
        database.execute("CREATE INDEX hit_rank ON hits(q,species,raw DESC,t)")
        count = database.execute("SELECT COUNT(*) FROM hits").fetchone()[0]
        if count <= 0:
            raise ValueError("No numeric hits")
        variants = {label: write_variant(database, output / label, names, species, cap)
                    for label, cap in (("all_hits", None), ("top100", 100))}
    if variants["all_hits"]["audit"]["summary"]["hits"] != count:
        raise ValueError("All-hit checkpoint lost rows")
    return variants, count


def convert(root, admission_path, admission_sha, admission_job, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if not os.environ.get("SLURM_JOB_ID", "").isdigit() or os.environ.get("SLURM_CPUS_PER_TASK") != "2":
        raise ValueError("Require scheduled two-CPU conversion")
    accounting = subprocess.check_output(["sacct", "-j", str(admission_job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS,ReqMem"], text=True)
    scheduler = require_completed_job(accounting, admission_job)
    if (scheduler["NodeList"], scheduler["AllocCPUS"]) != ("bizon", "2") or scheduler["ReqMem"] not in ("64G", "64Gn"):
        raise ValueError("Unexpected search-admission allocation")
    admission = read_frozen(admission_path, admission_sha)
    auditor = root / "benchmarks/work/publication_qfo_sequence_search_admission_v1"
    executor = root / "benchmarks/work/publication_qfo_sequence_search_v1"
    for path, revision in ((auditor, ADMITTER), (executor, EXECUTOR)):
        if subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip() != revision:
            raise ValueError("Frozen search source changed")
        subprocess.run(["git", "-C", str(path), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    if (admission["status"] != "corrected_qfo_search_panel_admitted_pending_numeric_validation"
            or admission["source"] != record(auditor / "benchmark_tools/admit_qfo_sequence_search_control.py")
            or admission["numeric_validated"] is not False or admission["accuracy_evaluated"] is not False):
        raise ValueError("Require unconverted corrected search admission")
    manifest = root / "benchmarks/work/qfo_sequence_search_control_v1/manifest.json"
    plan = verify_plan(manifest)
    execution_record = admission["execution"]
    check(execution_record)
    execution = read_frozen(Path(execution_record["path"]), execution_record["sha256"])
    records = validate_execution(execution, plan, admission["scheduler"], record(manifest), executor)
    if admission["targets"] != execution["targets"] or admission["manifest"] != record(manifest):
        raise ValueError("Admitted panel differs from execution")
    checked = [record(admission_path), execution_record, *records, *admission["checked_records"],
               plan["gene_metadata"], *plan["inputs"]]
    for item in checked:
        check(item)
    metadata = read_frozen(Path(plan["gene_metadata"]["path"]), plan["gene_metadata"]["sha256"])
    validate_metadata(plan["inputs"], metadata)
    output.mkdir(parents=True, exist_ok=False)
    report = {"status": "converting", "source": record(__file__), "admission": record(admission_path),
        "admission_scheduler": scheduler, "job_id": os.environ["SLURM_JOB_ID"], "checked_records": checked,
        "normalization": plan["normalization"], "accuracy_evaluated": False, "publication_ready": False,
        "sqlite_version": sqlite3.sqlite_version, "numpy_version": np.__version__,
        "helpers": [record(Path(__file__).with_name(n)) for n in (
            "convert_sequence_search_control.py", "audit_accuracy_checkpoint.py", "admit_qfo_sequence_search_control.py")],
        "checkpoint_writer": record(Path(__file__).resolve().parent.parent / "orthohmm/accuracy.py")}
    try:
        variants, count = convert_hits(plan, metadata, output)
        for item in [*checked, report["source"], *report["helpers"], report["checkpoint_writer"]]:
            check(item)
        report.update(status="corrected_qfo_sequence_numeric_checkpoints_verified", variants=variants, hits=count,
            genes=len(metadata), proteomes=78, numeric_validated=True,
            limitations=["Numeric identity/score integrity, not matched sensitivity or biological accuracy.",
                "Self hits and direction are retained; top100 is a post-search per-query/target-species diagnostic.",
                "Frozen graph replay and independent scoring remain required; no partial-panel inference."])
    except BaseException as error:
        report.update(status="failed", numeric_validated=False, error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "manifest.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("root", "admission", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    for flag in ("admission-sha256", "admission-job"):
        parser.add_argument("--" + flag, required=True)
    args = parser.parse_args()
    convert(args.root.resolve(), args.admission.resolve(), args.admission_sha256, args.admission_job, args.output.absolute())
