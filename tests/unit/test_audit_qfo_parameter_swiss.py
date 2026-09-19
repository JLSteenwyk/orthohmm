import gzip
from pathlib import Path

import pytest

from benchmark_tools import audit_qfo_parameter_swiss as module
from benchmark_tools.bootstrap_qfo_parameter_neighborhood import calculate
from tests.unit.test_audit_qfo_factorial_swiss import synthetic


def fixture(tmp_path):
    original, baseline = synthetic(tmp_path)
    entries = []
    for arm, entry in zip(module.ARMS, original):
        assessment = entry["assessment"]
        metrics = {r["metrics"]["metric_id"]: r["metrics"]["value"]
                   for r in assessment["native_assessments"] if r["challenge_id"] == "SwissTrees"}
        p, r = metrics["PPV"], metrics["TPR"]
        assessment["endpoints"] = {"SwissTrees": {"score": 2 * p * r / (p + r)}}
        entries.append({"arm": arm, "status": "admitted", "raw_file": entry["raw_file"], "assessment": assessment})
    return entries, baseline


def test_raw_counts_feed_paired_calculation(tmp_path):
    entries, baseline = fixture(tmp_path)
    counts = module.assemble(entries, baseline)
    assert counts["reference_relation_count"] == 270
    assert counts["arms"][0]["families"][0]["counts_without_prior"] == {"TP": 3, "FN": 5, "FP": 2, "TN": 5}
    result = calculate(counts, replicates=100)
    assert len(result["comparisons"]) == 6
    assert all(r["status"] == "estimated" for r in result["comparisons"])


def test_missing_arms_remain_explicit(tmp_path):
    entries, baseline = fixture(tmp_path)
    for index in (1, 2):
        entries[index] = {"arm": module.ARMS[index], "status": "not_admitted", "reason": "workflow not yet implemented"}
    counts = module.assemble(entries, baseline)
    assert counts["arms"][1:3] == entries[1:3]
    result = calculate(counts, replicates=100)
    assert result["multiplicity_endpoints"] == 18
    assert result["comparisons"][0]["metrics"] is None


@pytest.mark.parametrize("problem", ["order", "missing", "reference", "families", "coverage", "score", "truth", "members", "duplicate", "imputed"])
def test_changed_reference_or_metrics_rejected(tmp_path, problem):
    entries, baseline = fixture(tmp_path)
    if problem == "order":
        entries.reverse()
    elif problem == "missing":
        entries.pop()
    elif problem == "reference":
        baseline["reference"]["sha256"] = "other"
    elif problem == "families":
        entries[1]["assessment"]["swiss_reference_families"] = list(reversed(baseline["families"]))
    elif problem == "coverage":
        baseline["reference_orientation"]["family0"]["forward_relations"] += 1
    elif problem == "score":
        entries[1]["assessment"]["endpoints"]["SwissTrees"]["score"] += .01
    elif problem == "imputed":
        entries[1].update(status="not_admitted", reason="failed")
    else:
        path = Path(entries[1]["raw_file"]["path"])
        with gzip.open(path, "rt") as stream:
            lines = stream.readlines()
        if problem == "truth":
            a, b = lines[1].rstrip().split("\t"), lines[2].rstrip().split("\t")
            a[-1], b[-1] = b[-1], a[-1]
            lines[1], lines[2] = "\t".join(a) + "\n", "\t".join(b) + "\n"
        elif problem == "members":
            lines = [line.replace("family0_gene0", "changed_gene") for line in lines]
        else:
            lines.append(lines[1])
        with gzip.open(path, "wt") as stream:
            stream.writelines(lines)
    with pytest.raises(ValueError):
        module.assemble(entries, baseline)


def execution_fixture():
    report = {"scheduler": {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon",
              "AllocCPUS": "8", "JobIDRaw": "123"}, "index": 0, "pairs_manifest": {}, "conversion": {}}
    execution = {"status": "process_succeeded_pending_independent_admission", "exit_code": 0,
        "job_id": "123", "index": 0, "variant": "norm_low", "pairs_manifest": {}, "stage": {},
        "outputs": [{"path": "/results/SwissTrees/test.raw.txt.gz"}]}
    return report, execution


@pytest.mark.parametrize("problem", [None, "running", "exit", "job", "index", "variant", "pairs", "stage", "missing", "duplicate"])
def test_execution_binding(problem):
    report, execution = execution_fixture()
    if problem == "running":
        report["scheduler"]["State"] = "RUNNING"
    elif problem == "exit":
        execution["exit_code"] = True
    elif problem == "job":
        execution["job_id"] = "other"
    elif problem == "index":
        execution["index"] = 1
    elif problem == "variant":
        execution["variant"] = "norm_high"
    elif problem == "pairs":
        execution["pairs_manifest"] = {"path": "other"}
    elif problem == "stage":
        execution["stage"] = {"participant": "other"}
    elif problem == "missing":
        execution["outputs"] = []
    elif problem == "duplicate":
        execution["outputs"] *= 2
    if problem:
        with pytest.raises(ValueError):
            module.raw_from_execution("norm_low", report, execution)
    else:
        assert module.raw_from_execution("norm_low", report, execution) == execution["outputs"][0]


@pytest.mark.parametrize("arm", ["cpm_low", "cpm_high", "unknown"])
def test_unimplemented_admission_routes_fail_closed(arm):
    with pytest.raises(ValueError, match="workflow not yet"):
        module.validate_admission(arm, {}, {})
