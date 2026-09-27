"""Bound repository advisory ranges against the verified recovery installation."""

import argparse
import json
from pathlib import Path
import subprocess

from packaging.utils import canonicalize_name

from benchmark_tools.audit_dependency_lock import evaluate
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_publication_pipeline import save


def verify_inventory(installation, observed):
    expected = [(canonicalize_name(r["metadata"]["name"]), r["metadata"]["version"])
                for r in installation["install"]]
    actual = [(canonicalize_name(r["name"]), r["version"]) for r in observed]
    if len(set(expected)) != len(expected) or len(set(actual)) != len(actual) or sorted(expected) != sorted(actual):
        raise ValueError("Live distribution inventory differs from the installation report")
    return dict(package=[dict(name=n, version=v) for n, v in sorted(actual)])


def audit(repo, snapshot):
    directory = repo / "benchmarks/work/publication_recovery_install_20260926"
    receipt_path = repo / "benchmark_tools/results/publication_recovery_install_20260926_v2.json"
    receipt = json.loads(receipt_path.read_text())
    for key in ("install_report", "lock", "audit"):
        check(receipt[key])
    report_path = directory / "install_report.json"
    if record(report_path) != receipt["install_report"]:
        raise ValueError("Wrong recovery install report")
    watched = [record(snapshot), record(receipt_path), receipt["install_report"], receipt["lock"], receipt["audit"]]
    alerts = json.loads(snapshot.read_text())
    if (alerts["repository"] != "JLSteenwyk/orthohmm" or alerts["state_filter"] != "open"
            or any(a["state"] != "open" for a in alerts["alerts"])):
        raise ValueError("Wrong repository advisory snapshot scope")
    code = ("import json,importlib.metadata as m;"
            "print(json.dumps([dict(name=d.metadata['Name'],version=d.version) for d in m.distributions()]))")
    command = [str(directory / "venv/bin/python"), "-I", "-c", code]
    observed = json.loads(subprocess.check_output(command, text=True, timeout=60))
    inventory = verify_inventory(json.loads(report_path.read_text()), observed)
    comparisons = evaluate(inventory, alerts["alerts"])
    for item in watched:
        check(item)
    return dict(status="recovery_inventory_repository_advisory_ranges_checked", snapshot=record(snapshot),
        retrieved_utc=alerts["retrieved_utc"], inventory=inventory, command=command, comparisons=comparisons,
        affected_alerts=sum(r["status"] == "affected" for r in comparisons),
        alerts_without_installed_package=sum(r["status"] == "not_in_lock" for r in comparisons),
        checked_records=watched, source=record(__file__), evaluator=record(Path(evaluate.__code__.co_filename)),
        repository_alerts_closed=False, comprehensive_security_clearance=False, publication_ready=False,
        limitations=["Only ranges returned by the repository-alert snapshot, not all installed dependencies' vulnerabilities.",
                     "Live version inventory checked; prior byte audit is referenced, not rerun here.",
                     "No reachability, build-chain, OS/native-library or exploitability clearance.",
                     "Historical locks and repository alerts remain unchanged."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "snapshot", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.repo.resolve(), args.snapshot.resolve())
    save(args.output.resolve(), result)
    print(json.dumps({k: result[k] for k in ("status", "affected_alerts", "alerts_without_installed_package")}))
