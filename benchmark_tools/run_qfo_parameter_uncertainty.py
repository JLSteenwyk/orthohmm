"""Reconstruct admitted parameter-family counts and run the frozen bootstrap."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_parameter_swiss import audit
from benchmark_tools.bootstrap_qfo_parameter_neighborhood import calculate, REPLICATES, SEED, MULTIPLICITY
from benchmark_tools.cpm_replay_context import PLAN_SHA, PROTOCOL_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

SOURCES = {
    "audit_qfo_parameter_swiss.py": "dc166b360e8d94cdd3d29775bd19d178199f9c876fa1b288781624176fd2591a",
    "bootstrap_qfo_parameter_neighborhood.py": "369e5059dff46a6a97c0e01c28696f79b10c394bba1883d3e6a7878a329b149e",
    "audit_qfo_factorial_swiss.py": "50a2e722c0bd57fbf8fb96412ce3e4f788132491983563dcc75b26032ee13f28",
    "audit_qfo_swiss_counts.py": "ff8e5314c98bdffce54732df94611e9d340d854e70392358c495b636b67452ae",
    "bootstrap_qfo_swiss_stages.py": "0598773a8bd04846f24e4cb27ab7b3529410c11451cbb399b9148a320fd06435",
    "export_qfo_corrected_factorial.py": "8a508b229d6aee560a65c7700d4281f4b8fd550b748ef9062e18bc9f84717ecd",
    "run_qfo_parameter_assessment.py": "9da879481cf568c6037f29e2c726d24ce861fc99e55eaa55bed5b76e41b7cc3a",
    "run_qfo_cpm_assessment.py": "a8fad35977148e89557132f43362c9c842b602409de3974440dae93f6676c048",
    "validate_qfo_native_assessment.py": "4786081a21ce04b8ae25f1f3e83ed8a1766472ff7e4476ae70752fe7f0656da5",
    "prepare_ob_candidate_neighborhood.py": "03843ed17fa9c44ea1ce5249782cdf009d8c7db9ceb85c62aa85bc2b4156ec29",
    "run_simulation_methods.py": "cccf579c221c51e8d03d512264efc3d0a61f4c0a598863d2423bdbd90f7921be",
    "cpm_replay_context.py": "640da1894ae3b2e6a6b9db5ba8b060ef37c7fda2231a8591f57aa043a8b1e49b",
}


def run(inventory_path, inventory_sha, baseline_path, plan_path, protocol_path, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError("Require a fresh parameter uncertainty path")
    source = record(__file__)
    helpers = [record(Path(__file__).with_name(name)) for name in SOURCES]
    if any(item["sha256"] != SOURCES[name] for name, item in zip(SOURCES, helpers)):
        raise ValueError("Changed frozen parameter analysis implementation")
    protocol = record(protocol_path)
    if protocol["sha256"] != PROTOCOL_SHA:
        raise ValueError("Changed frozen parameter protocol")
    read_frozen(plan_path, PLAN_SHA)
    inventory_record = record(inventory_path)
    if inventory_record["sha256"] != inventory_sha:
        raise ValueError("Changed parameter admission inventory")
    checked = [source, *helpers, protocol, record(plan_path), inventory_record, record(baseline_path)]
    for item in checked:
        check(item)
    counts = audit(inventory_path, inventory_sha, baseline_path)
    if counts["source"] != record(Path(__file__).with_name("audit_qfo_parameter_swiss.py")):
        raise ValueError("Wrong reconstructed count-audit implementation")
    result = calculate(counts)
    if (result["protocol_controls_match"] is not True or result["replicates"] != REPLICATES
            or result["seed"] != SEED or result["multiplicity_endpoints"] != MULTIPLICITY):
        raise ValueError("Bootstrap controls differ from frozen protocol")
    checked.extend(counts["checked_inputs"])
    for item in checked:
        check(item)
    estimable = sum(row["status"] == "estimated" for row in result["comparisons"])
    result.update(status="corrected_qfo_parameter_uncertainty_audited", source=source, helpers=helpers,
        protocol=protocol, plan=record(plan_path), admission_inventory=inventory_record,
        baseline_audit=record(baseline_path), reconstructed_counts=counts, checked_inputs=checked,
        scientific_inputs_admitted=True, uncertainty_admitted=estimable > 0,
        estimated_contrasts=estimable, complete_panel=estimable == 6, publication_ready=False)
    result["limitations"] = [text for text in result["limitations"] if not text.startswith("Numerical calculation only;")]
    result["limitations"].extend([
        "Raw family predictions and native score arithmetic were rechecked; upstream inference/scoring admissions were not rerun.",
        "Uncertainty admission applies only to estimated contrasts; complete_panel is false if any planned contrast is unavailable.",
        "A control-only integration check produces no admitted variant intervals and does not complete parameter robustness."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("inventory", "baseline", "plan", "protocol", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--inventory-sha256", required=True)
    args = parser.parse_args()
    run(args.inventory.resolve(), args.inventory_sha256, args.baseline.resolve(), args.plan.resolve(),
        args.protocol.resolve(), args.output.absolute())
