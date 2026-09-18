from collections import Counter
from copy import deepcopy
import gzip
import itertools
import json
from pathlib import Path

import pytest

from benchmark_tools import audit_qfo_factorial_swiss as module
from benchmark_tools.audit_qfo_swiss_counts import HEADER, statistics, REFERENCE_SHA
from benchmark_tools.bootstrap_qfo_factorial import CELLS, validated_values
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.validate_qfo_native_assessment import AXES


def synthetic(tmp_path):
    families = [f"family{i}" for i in range(18)]
    entries = []
    for index, label in enumerate(CELLS):
        lines, scores = [], {}
        for family in families:
            counts = Counter()
            genes = [f"{family}_gene{i}" for i in range(6)]
            for pair_index, (a, b) in enumerate(itertools.combinations(genes, 2)):
                truth, predicted = pair_index % 2 == 0, (pair_index + index) % 3 == 0
                outcome = ("TP" if predicted else "FN") if truth else ("FP" if predicted else "TN")
                counts[outcome] += 1
                lines.append(f"{family}\t{a}\t{b}\t{outcome}\n")
            scores[family] = statistics(counts)
        raw = tmp_path / f"{label}.raw.txt.gz"
        with gzip.open(raw, "wt") as stream:
            stream.write(HEADER + "\n" + "".join(lines))
        participant = f"ohmm_qfo_factorial_{label}"
        metrics = []
        for challenge, axes in {**AXES, **{f"SwissTrees-{f}": ("PPV", "TPR") for f in families}}.items():
            for metric in axes:
                value = .5 if metric != "NR_ORTHOLOGS" else 10
                if challenge == "SwissTrees":
                    value = sum(scores[f][metric] for f in families) / 18
                elif challenge.startswith("SwissTrees-"):
                    value = scores[challenge.removeprefix("SwissTrees-")][metric]
                metrics.append({"type": "assessment", "community_id": "QfO", "participant_id": participant,
                    "_id": f"{challenge}-{metric}", "challenge_id": challenge,
                    "metrics": {"metric_id": metric, "value": value, "stderr": 0}})
        entries.append({"cell": label, "raw_file": record(raw), "assessment": {
            "participant": participant, "swiss_reference_families": families, "native_assessments": metrics}})
    baseline = {"status": "raw_swiss_family_counts_verified", "reference": {"sha256": REFERENCE_SHA},
                "families": families, "stages": [{"raw_file": entries[0]["raw_file"]}],
                "shared_represented_genes": {}, "reference_orientation": {
                    f: {"forward_relations": 15, "mapped_proteins": 6} for f in families}}
    return entries, baseline


def test_counts_reproduce_all_native_scores_and_bootstrap_input(tmp_path):
    entries, baseline = synthetic(tmp_path)
    result = module.assemble(entries, baseline)
    assert result["reference_relation_count"] == 270
    assert validated_values(result).shape == (8, 18, 2)
    assert result["cells"][0]["families"][0]["counts_without_prior"] == {"TP": 3, "FN": 5, "FP": 2, "TN": 5}


@pytest.mark.parametrize("mutation", ["missing", "order", "family", "score", "truth", "members", "coverage", "duplicate"])
def test_inconsistent_evidence_rejected(tmp_path, mutation):
    entries, baseline = synthetic(tmp_path)
    if mutation == "missing":
        entries.pop()
    elif mutation == "order":
        entries.reverse()
    elif mutation == "family":
        entries[1]["assessment"]["swiss_reference_families"] = list(reversed(baseline["families"]))
    elif mutation == "score":
        entries[1]["assessment"]["native_assessments"][-1]["metrics"]["value"] = .999
    elif mutation == "coverage":
        baseline["reference_orientation"]["family0"]["forward_relations"] = 16
    else:
        path = Path(entries[1]["raw_file"]["path"])
        with gzip.open(path, "rt") as stream:
            lines = stream.readlines()
        if mutation == "truth":
            # Swap a positive/negative truth label without changing aggregate counts.
            a, b = lines[1].rstrip().split("\t"), lines[2].rstrip().split("\t")
            a[-1], b[-1] = b[-1], a[-1]
            lines[1], lines[2] = "\t".join(a) + "\n", "\t".join(b) + "\n"
        elif mutation == "members":
            lines = [line.replace("family0_gene0", "family0_newgene") for line in lines]
        else:
            lines.append(lines[1])
        with gzip.open(path, "wt") as stream:
            stream.writelines(lines)
    with pytest.raises(ValueError):
        module.assemble(entries, baseline)


def baseline_report():
    return json.loads((Path(__file__).resolve().parents[2] /
        "benchmark_tools/results/qfo_factorial_assessment_p0_c0_r0_20260917.json").read_text())


def test_real_saved_baseline_binding():
    report = baseline_report()
    assert module.selected_admission(report, 0) == report["admitted_stage"]


@pytest.mark.parametrize("mutation", ["not_admitted", "index", "source", "pair", "participant"])
def test_reuse_binding_rejected(mutation):
    report = baseline_report()
    if mutation == "not_admitted":
        report["accuracy_admitted"] = False
    elif mutation == "index":
        report["index"] = 4
    elif mutation == "source":
        report["reused_admission"]["sha256"] = "0" * 64
    elif mutation == "pair":
        report["stage"]["pairs"]["sha256"] = "0" * 64
    else:
        report["original_participant"] = "wrong"
    with pytest.raises(ValueError):
        module.selected_admission(report, 0)


def test_fresh_admission_and_failure_rejection():
    report = {"cell": CELLS[1], "index": 1, "accuracy_admitted": True,
              "status": "fresh_factorial_assessment_admitted",
              "scheduler": {"State": "COMPLETED", "ExitCode": "0:0"},
              "assessment": {"participant": f"ohmm_qfo_factorial_{CELLS[1]}"}}
    assert module.selected_admission(report, 1) is report
    for key, value in (("status", "running"), ("scheduler", {"State": "FAILED", "ExitCode": "1:0"}),
                       ("assessment", {"participant": "wrong"})):
        changed = deepcopy(report)
        changed[key] = value
        with pytest.raises(ValueError):
            module.selected_admission(changed, 1)


def test_nested_file_inventory():
    item = {"path": "/tmp/example", "bytes": 7, "sha256": "abc"}
    assert list(module.file_records({"a": [item, {"b": item}], "c": "text"})) == [item, item]


def test_complete_file_audit_and_post_admission_tamper(tmp_path, monkeypatch):
    entries, baseline = synthetic(tmp_path)

    def save(name, value):
        path = tmp_path / name
        path.write_text(json.dumps(value))
        return record(path)

    native_rows = []
    conversions = []
    for index, entry in enumerate(entries):
        directory = tmp_path / str(index) / "SwissTrees"
        directory.mkdir(parents=True)
        raw = directory / "participant.raw.txt.gz"
        raw.write_bytes(Path(entry["raw_file"]["path"]).read_bytes())
        raw_record = record(raw)
        execution = save(f"execution{index}.json", {"status": "process_succeeded_pending_independent_admission",
            "exit_code": 0, "outputs": [raw_record]})
        conversion = {"stage": f"stage{index}", "pairs": raw_record, "filtered_pairs": raw_record,
            "partition": raw_record, "total_pairs": 10, "retained_pairs": 10, "removed_mapping_pairs": 0}
        conversions.append({**conversion, "candidate_partition": raw_record,
                            "reused_conversion": {"stage": conversion["stage"]}})
        native_rows.append({"status": "admitted", "index": 1 if index == 0 else 3,
            "participant": entry["assessment"]["participant"], "assessment": entry["assessment"],
            "conversion": conversion, "execution_report": execution})
    original = save("original.json", {"records": [None, native_rows[0], None, native_rows[4]]})
    monkeypatch.setattr(module, "ADMITTED_SHA", original["sha256"])
    manifest = {"cells": []}
    for index, native in enumerate(native_rows):
        report = {"cell": CELLS[index], "index": index, "accuracy_admitted": True}
        if index in (0, 4):
            report.update(status="admitted_reused_assessment", scoring_rerun=False,
                reused_admission=original, admitted_stage=native, stage=conversions[index],
                original_participant=native["participant"])
        else:
            report.update(status="fresh_factorial_assessment_admitted", assessment=native["assessment"],
                execution_report=native["execution_report"], scheduler={"State": "COMPLETED", "ExitCode": "0:0"})
        manifest["cells"].append({"cell": CELLS[index], "admission": save(f"admission{index}.json", report)})
    baseline_record = save("baseline.json", baseline)
    monkeypatch.setattr(module, "BASE_COUNTS_SHA", baseline_record["sha256"])
    inventory = save("inventory.json", manifest)
    result = module.audit(Path(inventory["path"]), inventory["sha256"], Path(baseline_record["path"]))
    assert len(result["cells"]) == 8
    assert len(result["checked_inputs"]) > 25
    with pytest.raises(ValueError, match="Changed admission inventory"):
        module.audit(Path(inventory["path"]), "0" * 64, Path(baseline_record["path"]))
    raw.write_bytes(b"changed after admission")
    with pytest.raises(ValueError, match="Frozen input/source identity changed"):
        module.audit(Path(inventory["path"]), inventory["sha256"], Path(baseline_record["path"]))
