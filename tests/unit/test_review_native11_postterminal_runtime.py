import copy
from datetime import datetime, timezone
import json

import pytest

from benchmark_tools import review_native11_postterminal_runtime as current
from benchmark_tools.inspect_native_python_lookup import compare_lookup
from benchmark_tools.prepare_allocated_native_factorial_qfo_pairs import admit_conversion


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))
    return current.record(path)


@pytest.fixture
def context(tmp_path):
    roots = [tmp_path / name for name in ("os", "private")]
    for root in roots:
        root.mkdir()
        (root / "original").write_text("original")
    manifests = [write(tmp_path / f"manifest_{i}.json", current.inventory([root]))
                 for i, root in enumerate(roots)]
    specs = [(r["path"], r["sha256"]) for r in manifests]
    expected = [dict(path=r["path"], sha256=r["sha256"], status="runtime_tree_identity_matches",
        records=2, scientific_execution_authorized=False) for r in manifests]
    baseline = write(tmp_path / "baseline.json", {})
    report = dict(python="3.10", executable="/fixture/python", cwd="/fixture",
        requested=["fixture"], paths=[], modules={}, mapped_files=[], meta_path=[],
        path_hooks=[], editable={}, dont_write_bytecode=True, files=[],
        coverage=dict(missing=[], changed=[], all_covered=True), scientific_origin={})
    reference = write(tmp_path / "expected_lookup.json", report)
    binding = write(tmp_path / "binding.json", dict(runtime_specs=specs))
    lookup = write(tmp_path / "lookup.json", dict(binding=binding, baseline=baseline,
        source=current.record(current.ROOT / "benchmark_tools/inspect_native_python_lookup.py"),
        interpreters={name: dict(reports=[reference]) for name in ("orthohmm", "orthofinder")}))
    session = tmp_path / "session"
    output = tmp_path / "output"
    original = tmp_path / "original.fa"
    original.write_text(">gene\nAAAA\n")
    copy_path = output / "input" / original.name
    copy_path.parent.mkdir(parents=True)
    copy_path.write_bytes(original.read_bytes())
    run = dict(inputs=[current.record(original)], native_order=[original.name], output_root=str(output))
    verification, phases = {}, {}
    for number, side in enumerate(("before", "after"), 1):
        comparisons = {}
        for name in ("orthohmm", "orthofinder"):
            observed = write(session / f"lookup_checks/check_{number:02d}/{name}.json", report)
            comparisons[name] = dict(compare_lookup(report, report), report=observed)
        checked = dict(status="runtime_and_lookup_checked", scientific_execution_authorized=False,
            runtime=expected, lookup=comparisons)
        write(session / f"lookup_checks/checked_{number:02d}.json", checked)
        prepared = None if side == "before" else dict(datasets=[dict(native_order=run["native_order"],
            inputs_in_native_order=[current.record(copy_path)])])
        verification[side] = dict(runtime=checked, original_inputs=run["inputs"], prepared_inputs=prepared)
        verification[side + "_check_wall_s"] = 1.0
        phases[side] = dict(lookup=comparisons, inventories=expected)
    prior = dict(status="runtime_brackets_and_lookup_replayed", phases=phases,
        fresh_terminal_inventory=expected, continuous_runtime_integrity_established=False)
    additions = []
    for name in ("late_a", "late_b"):
        path = roots[0] / name
        path.write_text(name)
        additions.append(dict(current.record(path), kind="file", mode=path.stat().st_mode & 0o7777))
    times = tuple(datetime(2026, 10, 8, hour, tzinfo=timezone.utc) for hour in (8, 9, 10))
    return dict(plan=dict(runtime_lookup=lookup, baseline=baseline), run=run,
        session=session, verification=verification, prior=prior,
        evidence=current.Evidence(), additions=additions, times=times, roots=roots,
        specs=specs, manifest_refs=manifests)


def replay(context):
    return current.replay_runtime(**{k: context[k] for k in
        ("plan", "run", "session", "verification", "prior", "evidence", "additions", "times")})


def test_original_bracket_kernel_and_fresh_inventory_are_replayed(context):
    report = replay(context)
    assert report["schema"] == current.SCHEMA
    assert report["historical_phases"] == context["prior"]["phases"]
    assert report["historical_first_terminal_inventory"] == context["prior"]["fresh_terminal_inventory"]
    assert report["current_original_inventory_equality"] is False
    assert report["current_private_inventory_equality"] is True
    assert report["current_original_entries_unchanged"] is True
    assert report["current_inventories"][0]["comparison"]["added"] == context["additions"]
    assert report["current_inventory_check_wall_s"] > 0
    assert "fresh_terminal_inventory" not in report
    assert "fresh_terminal_check_wall_s" not in report
    for key in ("full_review_admitted", "terminal_reviewed", "next_identity_authorized",
                "accuracy_evaluated", "resource_measurements_admitted", "automatic_retry",
                "publication_ready", "continuous_runtime_integrity_established",
                "original_ordinary_full_review_success"):
        assert report[key] is False
    assert len(context["evidence"].finish()) >= 12


@pytest.mark.parametrize("change", ["change_original", "remove_original", "unknown_addition",
                                  "change_addition", "private_change", "private_addition"])
def test_fresh_inventory_differences_are_not_hidden(context, change):
    os_root, private = context["roots"]
    if change == "change_original":
        (os_root / "original").write_text("changed")
    elif change == "remove_original":
        (os_root / "original").unlink()
    elif change == "unknown_addition":
        (os_root / "unexpected").write_text("extra")
    elif change == "change_addition":
        (os_root / "late_a").write_text("changed")
    elif change == "private_change":
        (private / "original").write_text("changed")
    else:
        (private / "extra").write_text("extra")
    with pytest.raises(ValueError):
        replay(context)


@pytest.mark.parametrize("change", ["before", "after", "lookup", "prior_phase", "prior_descriptor",
                                  "manifest", "copied_input", "chronology"])
def test_original_bracket_and_historical_binding_checks_remain_enforced(context, change):
    if change in ("before", "after"):
        context["verification"][change]["original_inputs"] = []
    elif change == "lookup":
        path = context["session"] / "lookup_checks/check_01/orthohmm.json"
        value = json.loads(path.read_text())
        value["modules"] = {"untrusted": "/untrusted"}
        ref = write(path, value)
        checked_path = context["session"] / "lookup_checks/checked_01.json"
        checked = json.loads(checked_path.read_text())
        checked["lookup"]["orthohmm"]["report"] = ref
        write(checked_path, checked)
        context["verification"]["before"]["runtime"] = checked
    elif change == "prior_phase":
        context["prior"] = copy.deepcopy(context["prior"])
        context["prior"]["phases"]["before"]["lookup"] = {}
    elif change == "prior_descriptor":
        context["prior"] = copy.deepcopy(context["prior"])
        context["prior"]["fresh_terminal_inventory"][0]["records"] = 99
    elif change == "manifest":
        write(current.Path(context["manifest_refs"][0]["path"]), {})
    elif change == "copied_input":
        (current.Path(context["run"]["output_root"]) / "input/original.fa").write_text("changed")
    else:
        context["times"] = tuple(reversed(context["times"]))
    with pytest.raises(ValueError):
        replay(context)


def test_metadata_and_duplicate_paths_are_rejected(context, monkeypatch):
    original_inventory = current.inventory
    def changed(roots):
        snapshot = original_inventory(roots)
        snapshot["exclusions"] = []
        return snapshot
    monkeypatch.setattr(current, "inventory", changed)
    with pytest.raises(ValueError, match="metadata"):
        replay(context)
    snapshot = original_inventory([context["roots"][0]])
    snapshot["records"].append(snapshot["records"][0])
    with pytest.raises(ValueError, match="Duplicate"):
        current.compare_inventory(snapshot, snapshot)


def chronology_fixture():
    failure = dict(schema="native11_full_review_runtime_refusal_observation_v1",
        native_job_id=23985, original_failed_review_job_id=23986, failed_review_job_id=24033,
        original_reviewer_error="Runtime inventory changed",
        accounting=dict(stdout=("23985|COMPLETED|0:0|2026-10-07T15:27:53|2026-10-08T08:56:00|"
            "17:28:07|64|128G|bizon|gpu|1560|native\n"
            "23986|FAILED|0:11|2026-10-08T08:56:00|2026-10-08T09:41:50|00:45:50|2|128G|bizon|gpu|360|review\n"
            "24033|FAILED|1:0|2026-10-08T10:50:42|2026-10-08T10:51:51|00:01:09|2|128G|bizon|gpu|360|review\n")),
        package_log_matching_lines={"/var/log/dpkg.log": [
            "2026-10-08 10:40:45 install lftp:amd64 <none> 4.9.2-2ubuntu1.1"]})
    classification = dict(prior_inventory_time_is_review_terminal_upper_bound=True,
        timezone="America/New_York", package_install_matching_line=failure["package_log_matching_lines"]["/var/log/dpkg.log"][0],
        result=dict(native_end="2026-10-08T08:56:00-04:00",
            prior_inventory_latest_at="2026-10-08T09:41:50-04:00", installed_at="2026-10-08T10:40:45-04:00"))
    return classification, failure


def test_chronology_replays_retained_accounting_and_package_line():
    classification, failure = chronology_fixture()
    times = current.chronology(classification, failure)
    assert times[0] < times[1] < times[2]
    assert all(t.utcoffset().total_seconds() == -14400 for t in times)


@pytest.mark.parametrize("change", ["time", "zone", "upper_bound", "duplicate_job", "outcome", "package"])
def test_chronology_rejects_unsupported_history(change):
    classification, failure = chronology_fixture()
    if change == "time":
        classification["result"]["native_end"] = "2026-10-08T08:55:00-04:00"
    elif change == "zone":
        classification["timezone"] = "UTC"
    elif change == "upper_bound":
        classification["prior_inventory_time_is_review_terminal_upper_bound"] = False
    elif change == "duplicate_job":
        failure["accounting"]["stdout"] += failure["accounting"]["stdout"].splitlines()[0] + "\n"
    elif change == "outcome":
        failure["accounting"]["stdout"] = failure["accounting"]["stdout"].replace("COMPLETED", "FAILED")
    else:
        failure["package_log_matching_lines"]["/var/log/dpkg.log"] = []
    with pytest.raises(ValueError):
        current.chronology(classification, failure)


def test_ordinary_conversion_gate_still_rejects_new_component(context):
    report = replay(context)
    run = dict(index=11, dataset="qfo_corrected", cell="p1_c1_r0", repeat=0)
    with pytest.raises(ValueError, match="Require bound successful"):
        admit_conversion(report, {}, dict(plan={}, amendment={}, job_id=23985), run)


def test_existing_destination_refuses_before_any_runtime_work(tmp_path, monkeypatch):
    destination = tmp_path / "exists"
    destination.mkdir()
    monkeypatch.setattr(current, "DEFAULT_DESTINATION", destination)
    monkeypatch.setattr(current, "record", lambda *args: pytest.fail("Should refuse before reading sources"))
    with pytest.raises(ValueError, match="fresh fixed"):
        current.review(destination, "untrusted")
