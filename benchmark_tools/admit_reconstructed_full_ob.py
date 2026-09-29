"""Independently read back reconstructed-base job 22376 against admitted 22337."""

import argparse
import json
from pathlib import Path
import subprocess

from benchmark_tools.admit_integrated_full_ob import read
from benchmark_tools.audit_frozen_overlay_install import local_install_wheels
from benchmark_tools.audit_installed_orthobench import read_root_hogs
from benchmark_tools.audit_recovery_advisories import verify_inventory
from benchmark_tools.audit_recovery_install import installed_payload
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.readback_canonical_ob_phylogeny import comparison, SCORER_SHA
from benchmark_tools.run_installed_orthobench import fasta_ids
from benchmark_tools.run_publication_pipeline import save
from benchmark_tools.score_orthobench_partition import score_partition
from benchmark_tools.validate_reader_upgrade import verify_lock
from benchmark_tools.verify_reconstructed_full_ob_execution import verify, RESULTS, JOB

BASELINE_SHA = "5d3a7a1d0b84926ae7872b161520d5323514d3c22316ed067b2fe294a4c12b15"
NATIVE_FILES = (
    "orthohmm_pairwise_orthologs.tsv", "orthohmm_pairwise_orthologs_confidence.tsv",
    "orthohmm_reconciliation_nodes.tsv", "orthohmm_hierarchical_orthogroups.tsv",
)


def report_pins(scientific):
    reports = scientific["reports"]
    if set(reports) != {"structure", "sequences", "events", "hierarchy"}:
        raise ValueError("Expected all four scientific readers")
    pins = {}
    for item in reports.values():
        check(item)
        for pin in read(Path(item["path"]))["checked_records"]:
            previous = pins.setdefault(pin["path"], pin)
            if previous != pin:
                raise ValueError("Conflicting readback file pins")
    return pins


def baseline():
    path = RESULTS / "integrated_full_ob_result_22337.json"
    receipt = record(path)
    if receipt["sha256"] != BASELINE_SHA:
        raise ValueError("Changed admitted baseline")
    value = read(path)
    if value["job_id"] != 22337 or value["status"] != "full_integrated_orthobench_independently_admitted":
        raise ValueError("Wrong baseline job or admission")
    for item in value["evidence"]:
        check(item)
    scientific = next(r for r in value["evidence"] if r["path"].endswith("/scientific/result.json"))
    pins = report_pins(read(Path(scientific["path"])))
    check(value["partitions"]["current"])
    return value, pins, receipt


def native_comparison(old, new, pins):
    results = {}
    for name in NATIVE_FILES:
        historical = record(old / name)
        if pins.get(historical["path"]) != historical:
            raise ValueError("Historical output is not bound to admitted readback: " + name)
        current = record(new / name)
        results[name] = dict(historical=historical, current=current,
            byte_equal=(historical["bytes"], historical["sha256"]) == (current["bytes"], current["sha256"]))
    return results


def package_audit(root, plan):
    command = plan["command"]
    value = lambda key: Path(command[command.index(key) + 1])
    assets = value("--assets")
    result = {}
    for name, wheels, lock in (
        ("inference", assets / "wheels", assets / "benchmark_tools/results/publication_recovery_requirements_20260926.txt"),
        ("reader", value("--reader-wheels"), value("--reader-lock")),
    ):
        report = read(root / (name + "_install.json"))
        rows = local_install_wheels(report, wheels)
        verify_lock(lock.read_text(), rows)
        probe = "import json,importlib.metadata as m;print(json.dumps([dict(name=d.metadata['Name'],version=d.version) for d in m.distributions()]))"
        observed = json.loads(subprocess.check_output(
            [str(root / name / "bin/python"), "-I", "-c", probe], text=True, timeout=60))
        result[name] = dict(inventory=verify_inventory(report, observed), packages=[
            dict(name=w["name"], **installed_payload(Path(w["wheel"]["path"]),
                 root / name / "lib/python3.10/site-packages")) for w in rows])
    return result


def audit(directory, output):
    if output.exists():
        raise FileExistsError(output)
    admitted = verify(directory)
    historical, pins, receipt = baseline()
    plan = read(directory / "plan.json")
    root = directory / "run"
    output.mkdir(parents=True)
    save(output / "execution.json", admitted)
    save(output / "package_audits.json", package_audit(root, plan))
    readers = plan["command"][plan["command"].index("--readers") + 1]
    code = ("import sys;from pathlib import Path;sys.path.insert(0," + repr(readers) + ");"
            "from benchmark_tools.audit_publication_pipeline import audit;"
            "audit(Path(" + repr(str(root / "native")) + "),Path(" + repr(str(output / "scientific")) + "))")
    environment = dict(HOME="/tmp", PATH="/usr/bin:/bin", LANG="C.UTF-8", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    with (output / "scientific.log").open("x") as log:
        subprocess.run([str(root / "reader/bin/python"), "-I", "-B", "-c", code],
            env=environment, cwd=output, stdout=log, stderr=subprocess.STDOUT, timeout=3600, check=True)
    independent = read(output / "scientific/result.json")
    recorded = read(root / "score.json")
    if independent["summary"] != recorded["summary"] or recorded["scientific_readback"] != record(root / "readback/result.json"):
        raise ValueError("Independent scientific readback differs from recorded result")
    data = read(root / "data.json")
    universe = fasta_ids([Path(r["path"]) for r in data["fasta"]])
    if len(universe) != 251378 or len(data["references"]) != 70:
        raise ValueError("Wrong full dataset dimensions")
    old = Path(historical["partitions"]["current"]["path"])
    new = root / "native/inference/orthohmm_phylogeny/orthohmm_root_hogs.tsv"
    paths = dict(historical=old, current=new)
    partition_records = {k: record(p) for k, p in paths.items()}
    partitions = {k: read_root_hogs(p, universe) for k, p in paths.items()}
    if (recorded["dataset"] != "orthobench" or recorded["genes"] != len(universe)
            or recorded["groups"] != len(partitions["current"])):
        raise ValueError("Recorded score dimensions differ")
    refs, uncertain = [{Path(r["path"]).name: set(Path(r["path"]).read_text().splitlines())
                       for r in data[role]} for role in ("references", "uncertain")]
    scorer = record(Path(score_partition.__code__.co_filename))
    if scorer["sha256"] != SCORER_SHA:
        raise ValueError("Changed frozen scorer")
    scores = {k: score_partition(v, refs, uncertain) for k, v in partitions.items()}
    if scores["current"] != recorded["score"] or scores["historical"] != historical["scores"]["current"]:
        raise ValueError("Recomputed score differs from its own recorded score")
    native = native_comparison(old.parent, new.parent, pins)
    current_pins = report_pins(independent)
    for item in [partition_records["current"], *(r["current"] for r in native.values())]:
        if current_pins.get(item["path"]) != item:
            raise ValueError("Current output differs from independent scientific readback")
    compared = comparison(partitions["historical"], partitions["current"], scores["historical"], scores["current"])
    if verify(directory) != admitted or baseline()[2] != receipt:
        raise ValueError("Execution or baseline changed during verification")
    for item in native.values():
        check(item["historical"])
        check(item["current"])
    for item in partition_records.values():
        check(item)
    result = dict(status="reconstructed_full_orthobench_independently_admitted", job_id=JOB,
        baseline=receipt, source=record(__file__), scorer=scorer, scores=scores,
        comparison=compared, native_outputs=native, summary=independent["summary"],
        partitions=partition_records,
        evidence=[record(output / n) for n in ("execution.json", "package_audits.json", "scientific/result.json")],
        historical_scores_replaced=False, controlled_timing=False, publication_ready=False,
        limitations=["Independent admission validates this execution, not equality to the baseline.",
            "Native comparisons are byte-level; mismatches require interpretation, not automatic retries.",
            "Root comparison is label invariant. No biological generalization or controlled timing claim.",
            "Wheel audits exclude generated metadata/bytecode and relocated non-site data.",
            "Not a complete OS dependency, cross-host restoration or redistribution-rights audit."])
    save(output / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(audit(args.directory.resolve(), args.output.resolve())["status"])
