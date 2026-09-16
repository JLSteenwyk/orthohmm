import copy

import pytest

from benchmark_tools.admit_qfo_affinity import IDENTITY_KEYS, affinity_arms, validate_panel


def panel():
    cpus = list(range(32))
    arms = affinity_arms(cpus)
    rows = []
    for index, (label, selected) in enumerate(arms):
        worker = {key: {} for key in IDENTITY_KEYS}
        worker.update(cpu_affinity=selected, requested_cpu_affinity=selected, inherited_cpu_affinity=cpus, inputs=[])
        rows.append({"index": index, "arm": label, "worker": worker, "execution": {"exit_code": 0},
                     "software_identity_excluding_affinity_equal": True, "software_identity_equal": index == 2})
    return {"status": "affinity_panel_complete", "accuracy_evaluated": False, "job_id": "21311",
            "affinity_panel": True, "parent_cpu_affinity": cpus, "graph_inputs": [], "repeats": rows,
            "planned_arms": [{"name": label, "cpu_affinity": selected} for label, selected in arms]}


def test_complete_panel_identity_admitted_without_assuming_equal_partitions():
    validate_panel(panel())


@pytest.mark.parametrize("change", ["partial", "job", "order", "affinity", "inherited", "library", "inputs", "exit", "flag"])
def test_changed_affinity_evidence_rejected(change):
    report = copy.deepcopy(panel())
    row = report["repeats"][2]
    if change == "partial":
        report["repeats"].pop()
    elif change == "job":
        report["job_id"] = "other"
    elif change == "order":
        row["arm"] = "all_cpus_1"
    elif change == "affinity":
        row["worker"]["cpu_affinity"] = [1]
    elif change == "inherited":
        row["worker"]["inherited_cpu_affinity"] = [0]
    elif change == "library":
        row["worker"]["native_libraries"] = {"changed": True}
    elif change == "inputs":
        row["worker"]["inputs"] = ["changed"]
    elif change == "exit":
        row["execution"]["exit_code"] = 1
    elif change == "flag":
        row["software_identity_equal"] = False
    with pytest.raises(ValueError):
        validate_panel(report)
