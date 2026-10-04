from copy import deepcopy
import hashlib
import json
from pathlib import Path

import pytest

from benchmark_tools import link_qfo_native_cost as cost


def execution_fixture():
    argv = ["/python", "-m", "orthohmm", "/inputs", "--refinement_profile", "default", "--stop", "infer"]
    native = {"status": "corrected_high_sensitivity_native_evidence_admitted", "scheduler": {
        "JobIDRaw": "21707", "State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "32"}}
    execution = {"job_id": "21707", "method": "orthohmm_high_sensitivity", "exit_code": 0,
        "status": "process_succeeded_pending_native_admission", "native_argv": argv, "cwd": "/core"}
    plan = {"methods": {"orthohmm_high_sensitivity": {"native_argv": argv, "search_reuse": False, "cwd": "/core"}}}
    values = {"wall_s": 10., "user_cpu_s": 200., "system_cpu_s": 1., "peak_process_tree_rss_bytes": 1024}
    metrics = dict(values, command=["/python", "/core/orthohmm/__main__.py", *argv[3:]], cwd="/core",
        status="complete", stages={k: dict(values) for k in cost.STAGES}, metadata=dict(cost.COMMON),
        counts={"species": 78, "genes": 984137, "orthogroups": 391908},
        rss_measurement="sampled_sum_of_linux_proc_tree_rss")
    return native, execution, plan, metrics, {"exit_status": 0}, argv


def test_execution_keeps_distinct_intervals_and_memory_scopes():
    args = execution_fixture()
    before = deepcopy(args)
    cost.validate_execution(*args)
    assert args == before


@pytest.mark.parametrize("change", [
    lambda a: a[0]["scheduler"].update(State="FAILED"),
    lambda a: a[0]["scheduler"].update(NodeList="dgx"),
    lambda a: a[0]["scheduler"].update(AllocCPUS="64"),
    lambda a: a[1].update(job_id="different"),
    lambda a: a[1].update(exit_code=1),
    lambda a: a[2]["methods"]["orthohmm_high_sensitivity"].update(search_reuse=True),
    lambda a: a[3].update(command=["different"]),
    lambda a: a[3]["metadata"].update(leiden_seed=0),
    lambda a: a[3]["counts"].update(genes=984138),
    lambda a: a[3]["stages"].update(phylogeny={}),
    lambda a: a[3].update(rss_measurement="lifetime_cgroup_peak"),
    lambda a: a[3].update(wall_s=float("nan")),
    lambda a: a[3].update(user_cpu_s=-1),
    lambda a: a[3].update(peak_process_tree_rss_bytes=1024.),
    lambda a: a[3]["stages"]["search"].update(wall_s=0),
    lambda a: a[4].update(exit_status=1),
    lambda a: a[5].extend(["--start", "search_res"]),
])
def test_execution_refuses_wrong_settings_reuse_or_accounting(change):
    args = execution_fixture()
    change(args)
    with pytest.raises(ValueError):
        cost.validate_execution(*args)


def join_fixture():
    inputs = [{"path": "/inputs/%d.fa" % i, "bytes": i, "sha256": str(i)} for i in range(78)]
    native = {"path": "/native.json", "bytes": 10, "sha256": "native"}
    primary = {"path": "/primary.json", "bytes": 10, "sha256": "primary"}
    replay_pin = {"path": "/replay.json", "bytes": 10, "sha256": "replay"}
    admission = {"path": "/admission.json", "bytes": 10, "sha256": "admission"}
    seed = {"path": "/seed.txt", "bytes": 10, "sha256": "partition"}
    candidate = dict(seed, path="/candidate.txt")
    prepared = {"status": "corrected_qfo_four_candidate_arms_prepared_unscored", "plan": replay_pin,
        "input_fastas": inputs, "admission": admission, "candidate_arms": {"p1_c0": {
            "candidate_expansion": False, "seed_partition": seed, "candidate_partition": candidate}}}
    replay = {"status": "corrected_checked_replay_admitted", "plan": replay_pin,
        "coverage": [{"label": "strict_profiles_refined", "output": seed}]}
    plan = {"admission": native, "primary_plan": primary, "input_fastas": inputs}
    execution = {"manifest": primary}
    selected = {"dataset": "QfO", "key": "orthohmm_high_sensitivity", "input_records": inputs,
        "qfo_stage_provenance": {"cell": cost.CELL, "candidate_partition": candidate, "replay_admission": admission}}
    costs = {"dataset": "Corrected QfO", "cell": cost.CELL, "full_pipeline_wall_s": None,
        "full_pipeline_cpu_s": None, "full_pipeline_peak_memory_bytes": None}
    return prepared, replay, plan, native, execution, selected, costs


def test_join_does_not_replace_original_cached_costs():
    args = join_fixture()
    before = deepcopy(args)
    cost.validate_join(*args)
    assert args == before


@pytest.mark.parametrize("change", [
    lambda a: a[0]["candidate_arms"]["p1_c0"].update(candidate_expansion=True),
    lambda a: a[0]["candidate_arms"]["p1_c0"]["candidate_partition"].update(sha256="changed"),
    lambda a: a[1].update(status="failed"),
    lambda a: a[2].update(primary_plan={}),
    lambda a: a[2].update(input_fastas=[]),
    lambda a: a[5]["qfo_stage_provenance"].update(cell="p1_c1_r1"),
    lambda a: a[5]["qfo_stage_provenance"].update(replay_admission={}),
    lambda a: a[6].update(full_pipeline_wall_s=100.),
])
def test_join_refuses_different_scientific_evidence_or_imputed_original_costs(change):
    args = join_fixture()
    change(args)
    with pytest.raises(ValueError):
        cost.validate_join(*args)


@pytest.mark.parametrize("records", [[], [{"path": "/other"}], [
    {"path": "/wanted", "sha256": "one"}, {"path": "/wanted", "sha256": "two"}]])
def test_missing_or_conflicting_admitted_pin_is_refused(records):
    with pytest.raises(ValueError):
        cost.unique_pin(records, "/wanted")


def test_duplicate_identical_admission_references_do_not_create_repeats():
    pin = {"path": "/wanted", "bytes": 4, "sha256": "one"}
    assert cost.unique_pin([pin, pin], "/wanted") == pin


def test_named_and_unnamed_partitions_compare_membership_not_labels(tmp_path):
    named, plain = tmp_path / "named.txt", tmp_path / "plain.txt"
    named.write_text("OG0: a b\nOG1: c\n")
    plain.write_text("c\nb a\n")
    assert cost.partition(named, "named_groups") == cost.partition(plain, "space_separated_groups")
    plain.write_text("a b\na c\n")
    with pytest.raises(ValueError, match="Duplicate gene"):
        cost.partition(plain, "space_separated_groups")


def test_existing_output_refused_before_collection(tmp_path, monkeypatch):
    monkeypatch.setattr("sys.argv", ["cost", "--root", str(tmp_path), "--output-directory", str(tmp_path)])
    monkeypatch.setattr(cost, "collect", lambda _: pytest.fail("Must not collect again"))
    with pytest.raises(ValueError, match="existing output"):
        cost.main()


def test_actual_association_against_raw_retained_inputs_partitions_and_metrics():
    root = Path(__file__).resolve().parents[2]
    report = json.loads((root / "benchmark_tools/results/qfo_native_configuration_cost_20261004/association.json").read_text())
    assert report["schema"] == "qfo_native_configuration_cost_v1"
    assert report["job_id"] == 21707 and report["cell"] == "p1_c0_r0" and report["repeats"] == 1
    assert len(report["checked_records"]) == 133 and len(report["input_fastas"]) == 78
    for pin in report["checked_records"]:
        data = Path(pin["path"]).read_bytes()
        assert len(data) == pin["bytes"]
        assert hashlib.sha256(data).hexdigest() == pin["sha256"]
    inputs = [line[1:].split()[0] for pin in report["input_fastas"]
              for line in Path(pin["path"]).read_text().splitlines() if line.startswith(">")]
    assert len(inputs) == len(set(inputs)) == 984137
    native_rows = [line.split(":", 1) for line in Path(report["partition"]["native"]["path"]).read_text().splitlines()]
    assert len({row[0] for row in native_rows}) == len(native_rows) == 391908
    native = {tuple(sorted(row[1].split())) for row in native_rows}
    factorial = {tuple(sorted(line.split())) for line in
                 Path(report["partition"]["factorial"]["path"]).read_text().splitlines()}
    assert len(native) == len(factorial) == 391908 and native == factorial
    assert sum(map(len, native)) == 984137
    assert {gene for group in native for gene in group} == set(inputs)
    metrics = json.loads(Path(report["native_metrics"]["path"]).read_text())
    assert report["metrics"] == {key: metrics[key] for key in report["metrics"]}
    assert report["recorded_stages"] == metrics["stages"]
    assert report["metrics"]["wall_s"] == 70886.239038
    assert report["metrics"]["peak_process_tree_rss_bytes"] == 18245820416
    companion = report["companion"]["measurement"]
    assert companion["elapsed_seconds"] == 70888. and companion["max_process_rss_kib"] == 11885516
    assert all(v is None for v in report["original_cached_full_costs"].values())
    original = json.loads(Path(report["sources"]["costs"]["path"]).read_text())
    assert len(original["rows"]) == 16
    assert all(row["full_pipeline_wall_s"] is None for row in original["rows"])
    assert report["native_inference_repeated"] is False and report["scoring_repeated"] is False
    assert report["controlled_comparative_resources"] is False and report["publication_ready"] is False
