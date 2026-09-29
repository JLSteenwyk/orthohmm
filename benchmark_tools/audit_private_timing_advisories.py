"""Check a private timing candidate against retained repository alert ranges."""

import argparse
import json
from pathlib import Path
import subprocess

from benchmark_tools.audit_dependency_lock import evaluate
from benchmark_tools.audit_recovery_advisories import verify_inventory
from benchmark_tools.build_private_timing_environment import record, write


def validate_snapshot(snapshot):
    alerts = snapshot["alerts"]
    if (snapshot["repository"] != "JLSteenwyk/orthohmm"
            or snapshot["state_filter"] != "open" or not alerts
            or any(a["state"] != "open" or a["dependency"]["package"]["ecosystem"] != "pip"
                   for a in alerts)
            or len({a["number"] for a in alerts}) != len(alerts)):
        raise ValueError("Require nonempty, unique, open pip alerts for this repository")
    return alerts


def compare(selected, observed, snapshot, previous):
    installation = {"install": [{"metadata": {"name": n, "version": v}}
                                for n, v in selected.items()]}
    inventory = verify_inventory(installation, observed)
    alerts = validate_snapshot(snapshot)
    old = validate_snapshot(previous)
    comparisons = evaluate(inventory, alerts)
    added = sorted({a["number"] for a in alerts} - {a["number"] for a in old})
    affected = [r for r in comparisons if r["status"] == "affected"]
    return dict(inventory=inventory, comparisons=comparisons, added_alert_numbers=added,
                affected_manifest_alerts=len(affected),
                affected_unique_advisories=len({(r["package"], r["ghsa_id"]) for r in affected}),
                unique_advisories=len({(r["package"], r["ghsa_id"]) for r in comparisons}))


def audit(python, candidate, snapshot, previous):
    watched = [record(p) for p in (python, candidate, snapshot, previous)]
    receipt = json.loads(candidate.read_text())
    if receipt["status"] not in ("private_packaging_historical_payload_aligned",
                                  "private_timing_environment_candidate_installed"):
        raise ValueError("Require an installed private candidate receipt")
    code = ("import json,importlib.metadata as m;"
            "print(json.dumps([dict(name=d.metadata['Name'],version=d.version) for d in m.distributions()]))")
    command = [str(python), "-I", "-B", "-c", code]
    observed = json.loads(subprocess.check_output(command, text=True, timeout=60))
    result = compare(receipt["selected"], observed, json.loads(snapshot.read_text()),
                     json.loads(previous.read_text()))
    if watched != [record(p) for p in (python, candidate, snapshot, previous)]:
        raise ValueError("Inputs changed during audit")
    result.update(command=command, inputs=watched, source=record(Path(__file__)),
                  evaluator=record(Path(evaluate.__code__.co_filename)),
                  inventory_verifier=record(Path(verify_inventory.__code__.co_filename)),
                  status="private_candidate_advisory_ranges_checked",
                  scientific_execution_authorized=False, comprehensive_security_clearance=False,
                  limitations=["Version metadata, not package payload or exploitability verification.",
                               "Only retained repository advisory ranges; missing manifest alerts do not imply safety.",
                               "No installation, alert dismissal, or modification of historical locks."])
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("python", "candidate", "snapshot", "previous", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.python.absolute(), args.candidate.resolve(), args.snapshot.resolve(),
                   args.previous.resolve())
    write(args.output, result)
    print(json.dumps({k: result[k] for k in ("status", "affected_manifest_alerts",
                                           "affected_unique_advisories", "added_alert_numbers")}))
