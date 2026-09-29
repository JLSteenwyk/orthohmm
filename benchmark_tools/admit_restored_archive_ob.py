"""Independently validate restored-archive job 22377 against admitted job 22376."""

import argparse
from pathlib import Path
import subprocess

from benchmark_tools.admit_integrated_full_ob import read
from benchmark_tools.admit_reconstructed_full_ob import (
    native_comparison, package_audit, report_pins, NATIVE_FILES,
)
from benchmark_tools.audit_installed_orthobench import read_root_hogs
from benchmark_tools.readback_canonical_ob_phylogeny import comparison, SCORER_SHA
from benchmark_tools.run_installed_orthobench import fasta_ids
from benchmark_tools.run_integrated_full_job import record, save
from benchmark_tools.run_integrated_publication_workflow import check
from benchmark_tools.score_orthobench_partition import score_partition
from benchmark_tools.verify_restored_archive_execution import verify, RESULTS, JOB

BASELINE_SHA = "81c4de02279aa08500c9275a611cf180a48f64c3891234848a1c6a0f81da0fc0"


def reproduction_equal(compared, native):
    if set(native) != set(NATIVE_FILES):
        raise ValueError("Require all four native comparisons")
    return (compared["partitions"]["label_invariant_equal"] and compared["score_objects_equal"]
            and all(r["byte_equal"] for r in native.values()))


def baseline():
    path = RESULTS / "reconstructed_full_ob_result_22376.json"
    pin = record(path)
    if pin["sha256"] != BASELINE_SHA:
        raise ValueError("Changed admitted baseline")
    value = read(path)
    if value["job_id"] != 22376 or value["status"] != "reconstructed_full_orthobench_independently_admitted":
        raise ValueError("Wrong admitted baseline")
    for item in value["evidence"]:
        check(item)
    scientific = next(p for p in value["evidence"] if p["path"].endswith("/scientific/result.json"))
    pins = report_pins(read(Path(scientific["path"])))
    check(value["partitions"]["current"])
    return value, pins, pin


def audit(directory, output):
    if output.exists():
        raise FileExistsError(output)
    execution = verify(directory)
    historical, pins, baseline_pin = baseline()
    plan = read(directory / "plan.json")
    root = directory / "run"
    output.mkdir(parents=True)
    save(output / "execution.json", execution)
    save(output / "package_audits.json", package_audit(root, plan))
    readers = plan["command"][plan["command"].index("--readers") + 1]
    code = ("import sys;from pathlib import Path;sys.path.insert(0," + repr(readers) + ");"
            "from benchmark_tools.audit_publication_pipeline import audit;"
            "audit(Path(" + repr(str(root / "native")) + "),Path(" + repr(str(output / "scientific")) + "))")
    env = dict(HOME="/tmp", PATH="/usr/bin:/bin", LANG="C.UTF-8", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    with (output / "scientific.log").open("x") as log:
        subprocess.run([str(root / "reader/bin/python"), "-I", "-B", "-c", code],
                       env=env, cwd=output, stdout=log, stderr=subprocess.STDOUT, timeout=3600, check=True)
    independent = read(output / "scientific/result.json")
    recorded = read(root / "score.json")
    if (independent["summary"] != recorded["summary"]
            or recorded["scientific_readback"] != record(root / "readback/result.json")):
        raise ValueError("Independent readback differs")
    data = read(root / "data.json")
    universe = fasta_ids([Path(p["path"]) for p in data["fasta"]])
    if len(universe) != 251378 or len(data["references"]) != 70:
        raise ValueError("Wrong full dataset dimensions")
    old = Path(historical["partitions"]["current"]["path"])
    new = root / "native/inference/orthohmm_phylogeny/orthohmm_root_hogs.tsv"
    paths = dict(historical=old, current=new)
    partition_pins = {k: record(p) for k, p in paths.items()}
    partitions = {k: read_root_hogs(p, universe) for k, p in paths.items()}
    if (recorded["dataset"] != "orthobench" or recorded["genes"] != len(universe)
            or recorded["groups"] != len(partitions["current"])):
        raise ValueError("Recorded dimensions differ")
    refs, uncertain = [{Path(p["path"]).name: set(Path(p["path"]).read_text().splitlines())
                       for p in data[role]} for role in ("references", "uncertain")]
    scorer = record(Path(score_partition.__code__.co_filename))
    if scorer["sha256"] != SCORER_SHA:
        raise ValueError("Changed frozen scorer")
    scores = {k: score_partition(v, refs, uncertain) for k, v in partitions.items()}
    if scores["current"] != recorded["score"] or scores["historical"] != historical["scores"]["current"]:
        raise ValueError("Recomputed scores differ from recorded scores")
    native = native_comparison(old.parent, new.parent, pins)
    current_pins = report_pins(independent)
    for pin in [partition_pins["current"], *(r["current"] for r in native.values())]:
        if current_pins.get(pin["path"]) != pin:
            raise ValueError("Current output not bound to independent readback")
    compared = comparison(partitions["historical"], partitions["current"], scores["historical"], scores["current"])
    if verify(directory) != execution or baseline()[2] != baseline_pin:
        raise ValueError("Execution or baseline changed during audit")
    for item in native.values():
        check(item["historical"])
        check(item["current"])
    for pin in partition_pins.values():
        check(pin)
    equivalent = reproduction_equal(compared, native)
    result = dict(status="restored_archive_orthobench_independently_audited", job_id=JOB,
        reproduction_equal=equivalent, baseline=baseline_pin, source=record(__file__), scorer=scorer,
        scores=scores, comparison=compared, native_outputs=native, partitions=partition_pins,
        summary=independent["summary"],
        evidence=[record(output / n) for n in ("execution.json", "package_audits.json", "scientific/result.json")],
        historical_scores_replaced=False, controlled_timing=False, publication_ready=False,
        limitations=["Execution validation and reproduction equality are separate; mismatches are retained.",
            "Same-host restoration, not new biological validation, controlled timing or OS closure.",
            "Wheel audit excludes generated metadata/bytecode and relocated non-site data.",
            "No redistribution-rights or comprehensive security clearance."])
    save(output / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(audit(args.directory.resolve(), args.output.resolve())["status"])
