from collections import Counter
from copy import deepcopy
import gzip
import itertools
import json
from pathlib import Path

import pytest

from benchmark_tools import audit_qfo_sequence_swiss as module
from benchmark_tools.audit_qfo_swiss_counts import HEADER, statistics
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.validate_qfo_native_assessment import AXES
from tests.unit.test_run_qfo_sequence_assessment import stage
from tests.unit.test_audit_qfo_corrected_factorial_swiss import execution_fixture
from tests.unit.test_export_qfo_corrected_factorial import fixture as hmm_report
from benchmark_tools import bootstrap_qfo_sequence


def sequence_reports():
    conversion = stage()
    scheduler = dict(State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="2", JobIDRaw="1")
    report = dict(status="corrected_sequence_assessment_admitted", variant="all_hits",
        accuracy_admitted=True, publication_ready=False, conversion=conversion, conversion_scheduler=scheduler,
        assessment={"participant": conversion["participant"]}, pairs_manifest={"path": "/pairs.json"},
        scheduler=dict(State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="8", JobIDRaw="2"))
    execution = dict(status="process_succeeded_pending_independent_admission", exit_code=0, job_id="2",
        variant="all_hits", stage=deepcopy(conversion), conversion_scheduler=deepcopy(scheduler),
        pairs_manifest=report["pairs_manifest"], outputs=[{"path": "/native/SwissTrees/test.raw.txt.gz"}])
    return report, execution, conversion


@pytest.mark.parametrize("problem", [None, "status", "variant", "accuracy", "publication", "running", "job",
                                    "stage", "conversion_job", "pairs", "missing", "duplicate", "boolean_exit"])
def test_sequence_execution_binding(problem):
    report, execution, conversion = sequence_reports()
    if problem in ("status", "variant"):
        report[problem] = "wrong"
    elif problem in ("accuracy", "publication"):
        report["accuracy_admitted" if problem == "accuracy" else "publication_ready"] = problem != "accuracy"
    elif problem == "running":
        report["scheduler"]["State"] = "RUNNING"
    elif problem == "job":
        execution["job_id"] = "3"
    elif problem == "stage":
        execution["stage"]["total_pairs"] = 9
    elif problem == "conversion_job":
        execution["conversion_scheduler"]["JobIDRaw"] = "3"
    elif problem == "pairs":
        execution["pairs_manifest"] = {}
    elif problem == "missing":
        execution["outputs"] = []
    elif problem == "duplicate":
        execution["outputs"] *= 2
    elif problem == "boolean_exit":
        execution["exit_code"] = False
    if problem:
        with pytest.raises(ValueError):
            module.selected_raw(report, execution, conversion, "all_hits")
    else:
        assert module.selected_raw(report, execution, conversion, "all_hits") == execution["outputs"][0]


def test_hmm_must_be_initial_only():
    report, execution = execution_fixture()
    with pytest.raises(ValueError, match="initial HMM"):
        module.selected_raw(report, execution, report["conversion"], "p0_c0_r0")


def synthetic(tmp_path):
    families = [f"family{i}" for i in range(18)]
    entries, orientation = [], {}
    for index, variant in enumerate(module.VARIANTS):
        lines, scores = [], {}
        for i, family in enumerate(families):
            genes = [f"{family}_gene{j}" for j in range(37 if i == 0 else 35)]
            pairs = list(itertools.combinations(genes, 2))[:650 if i == 0 else 595]
            orientation[family] = dict(forward_relations=len(pairs), mapped_proteins=len(genes))
            counts = Counter()
            for k, (a, b) in enumerate(pairs):
                truth, predicted = k % 2 == 0, (k + index) % 3 == 0
                outcome = ("TP" if predicted else "FN") if truth else ("FP" if predicted else "TN")
                counts[outcome] += 1
                lines.append(f"{family}\t{a}\t{b}\t{outcome}\n")
            scores[family] = statistics(counts)
        raw = tmp_path / variant / "SwissTrees/test.raw.txt.gz"
        raw.parent.mkdir(parents=True)
        with gzip.open(raw, "wt") as stream:
            stream.write(HEADER + "\n" + "".join(lines))
        metrics = []
        for challenge, axes in {**AXES, **{f"SwissTrees-{f}": ("PPV", "TPR") for f in families}}.items():
            for metric in axes:
                value = .5 if metric != "NR_ORTHOLOGS" else 10
                if challenge == "SwissTrees":
                    value = sum(scores[f][metric] for f in families) / 18
                elif challenge.startswith("SwissTrees-"):
                    value = scores[challenge.removeprefix("SwissTrees-")][metric]
                metrics.append(dict(type="assessment", community_id="QfO", participant_id=variant,
                    _id=f"{challenge}-{metric}", challenge_id=challenge,
                    metrics=dict(metric_id=metric, value=value, stderr=0)))
        p, r = (sum(scores[f][m] for f in families) / 18 for m in ("PPV", "TPR"))
        entries.append(dict(variant=variant, raw_file=record(raw), assessment=dict(participant=variant,
            swiss_reference_families=families, native_assessments=metrics,
            endpoints={"SwissTrees": {"score": 2*p*r/(p+r)}})))
    baseline = dict(status="raw_swiss_family_counts_verified", reference={"sha256": module.REFERENCE_SHA},
        families=families, stages=[{"raw_file": entries[0]["raw_file"]}],
        reference_orientation=orientation, shared_represented_genes={})
    return entries, baseline


def test_exact_full_reference_synthetic_assembly(tmp_path):
    entries, baseline = synthetic(tmp_path)
    result = module.assemble(entries, baseline)
    assert result["reference_relation_count"] == 10765
    assert module.validated_values(result).shape == (3, 18, 2)
    assert result["uncertainty_admitted"] is False
    assert result["variants"][0]["families"][0]["counts_without_prior"] != result["variants"][1]["families"][0]["counts_without_prior"]


@pytest.mark.parametrize("problem", ["missing", "order", "score", "reference", "coverage", "truth_swap", "duplicate"])
def test_invalid_raw_evidence(tmp_path, problem):
    entries, baseline = synthetic(tmp_path)
    if problem == "missing":
        entries.pop()
    elif problem == "order":
        entries.reverse()
    elif problem == "score":
        entries[1]["assessment"]["endpoints"]["SwissTrees"]["score"] += .01
    elif problem == "reference":
        baseline["reference"]["sha256"] = "changed"
    elif problem == "coverage":
        baseline["reference_orientation"]["family0"]["forward_relations"] += 1
    else:
        path = Path(entries[1]["raw_file"]["path"])
        with gzip.open(path, "rt") as stream:
            lines = stream.readlines()
        if problem == "truth_swap":
            a, b = (line.rstrip().split("\t") for line in lines[1:3])
            a[-1], b[-1] = b[-1], a[-1]
            lines[1:3] = ["\t".join(row) + "\n" for row in (a, b)]
        else:
            lines.append(lines[1])
        with gzip.open(path, "wt") as stream:
            stream.writelines(lines)
    with pytest.raises(ValueError):
        module.assemble(entries, baseline)


def test_source_bound_audit_and_changed_raw_file(tmp_path, monkeypatch):
    entries, baseline = synthetic(tmp_path)

    def save(name, value):
        path = tmp_path / name
        path.write_text(json.dumps(value))
        return record(path)

    marker = save("reference.json", {})
    baseline["reference"] = marker
    monkeypatch.setattr(module, "REFERENCE_SHA", marker["sha256"])
    monkeypatch.setattr(bootstrap_qfo_sequence, "REFERENCE_SHA", marker["sha256"])
    base = save("baseline.json", baseline)
    monkeypatch.setattr(module, "BASE_COUNTS_SHA", base["sha256"])
    inventory = {"variants": []}
    for index, entry in enumerate(entries):
        variant = entry["variant"]
        if index == 0:
            report = hmm_report(0)
            conversion = report["conversion"]
            conversion.update(pairs=marker, filtered_pairs=marker, native_admission_recheck=marker,
                              checked_records=[marker])
        else:
            report, _, conversion = sequence_reports()
            conversion.update(variant=variant, participant="ohmm_qfo_corrected_sequence_" + variant,
                pairs=marker, filtered_pairs=marker, graph_admission=marker, prediction=marker,
                checked_records=[marker])
            report["variant"] = variant
        participant = conversion["participant"]
        native = deepcopy(entry["assessment"])
        native["participant"] = participant
        for item in native["native_assessments"]:
            item["participant_id"] = participant
        if index == 0:
            report["assessment"].update({key: value for key, value in native.items() if key != "endpoints"})
            endpoint = report["assessment"]["endpoints"]["SwissTrees"]
            values = {m["metrics"]["metric_id"]: m["metrics"]["value"] for m in native["native_assessments"]
                      if m["challenge_id"] == "SwissTrees"}
            endpoint["native_participant"].update(metric_x=values["TPR"], metric_y=values["PPV"])
            endpoint["score"] = native["endpoints"]["SwissTrees"]["score"]
            report["assessment"]["secondary_six_metric_mean"] = (2.5 + endpoint["score"]) / 6
        else:
            report["assessment"] = native
        report["pairs_manifest"] = save(f"pairs{index}.json", conversion)
        report["scheduler"] = dict(State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="8", JobIDRaw="2")
        execution = dict(status="process_succeeded_pending_independent_admission", exit_code=0, job_id="2",
            stage=conversion, pairs_manifest=report["pairs_manifest"], conversion_scheduler=report["conversion_scheduler"],
            outputs=[entry["raw_file"]])
        execution.update({"index": 0, "cell": variant} if index == 0 else {"variant": variant})
        report["execution_report"] = save(f"execution{index}.json", execution)
        report["checked_records"] = [report["execution_report"]]
        inventory["variants"].append(dict(variant=variant, admission=save(f"admission{index}.json", report)))
    inv = save("inventory.json", inventory)
    result = module.audit(Path(inv["path"]), inv["sha256"], Path(base["path"]))
    assert result["reference_relation_count"] == 10765
    assert len(result["checked_inputs"]) >= 15
    assert result["source"] == record(module.__file__)
    Path(entries[2]["raw_file"]["path"]).write_bytes(b"changed after admission")
    with pytest.raises(ValueError, match="identity changed"):
        module.audit(Path(inv["path"]), inv["sha256"], Path(base["path"]))
