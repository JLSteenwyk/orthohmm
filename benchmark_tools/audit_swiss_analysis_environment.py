"""Compare the tested analysis environment with pins and retained advisory ranges."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from packaging.requirements import Requirement
from packaging.utils import canonicalize_name
from benchmark_tools.audit_dependency_lock import evaluate, identity


def check_pins(text, installed):
    expected = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        requirement = Requirement(line)
        specs = list(requirement.specifier)
        name = canonicalize_name(requirement.name)
        if (requirement.extras or requirement.marker or requirement.url or len(specs) != 1
                or specs[0].operator != "==" or "*" in specs[0].version or name in expected):
            raise ValueError("Require unique unconditional exact pins")
        expected[name] = specs[0].version
    actual = {canonicalize_name(name): version for name, version in installed.items()}
    if not expected or actual != expected:
        raise ValueError("Tested environment differs from exact pin inventory")
    return {"package": [{"name": name, "version": version} for name, version in sorted(actual.items())]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("requirements", "lock", "reproduction", "alerts", "output"):
        parser.add_argument("--" + flag, required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    evidence = [identity(p) for p in (args.requirements, args.lock, args.reproduction, args.alerts)]
    reproduced = json.loads(args.reproduction.read_text())
    if reproduced["status"] != "relocated_swiss_statistics_and_plot_workflow_reproduced":
        raise ValueError("Require successful reproduction evidence")
    packages = check_pins(args.requirements.read_text(), reproduced["environment"]["packages"])
    alerts = json.loads(args.alerts.read_text())["alerts"]
    if not alerts or any(r["dependency"]["manifest_path"] != "benchmark_tools/swiss_analysis_requirements.txt" for r in alerts):
        raise ValueError("Require a nonempty analysis-environment alert snapshot")
    rows = evaluate(packages, alerts)
    if any(r["status"] == "not_in_lock" for r in rows):
        raise ValueError("Advisory dependency missing from tested environment")
    if evidence != [identity(p) for p in (args.requirements, args.lock, args.reproduction, args.alerts)]:
        raise ValueError("Evidence changed during audit")
    result = {"status": "analysis_environment_ranges_evaluated", "evidence": evidence,
              "source": identity(Path(__file__)), "range_evaluator": identity(Path(__file__).with_name("audit_dependency_lock.py")),
              "comparisons": rows, "affected_alerts": sum(r["status"] == "affected" for r in rows),
              "limitations": ["Only retained advisory ranges checked, not a comprehensive security or exploitability audit.",
                              "Package versions are taken from the recorded successful reproduction environment.",
                              "Distribution hash enforcement is performed by uv pip sync, not reimplemented here.",
                              "Remote alert closure must be checked separately after push."]}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
