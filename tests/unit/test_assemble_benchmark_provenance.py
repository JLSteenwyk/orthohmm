import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import assemble_benchmark_provenance as provenance


def write_json(path, value):
    path.write_text(json.dumps(value))
    return provenance.record(path)


def execution(tmp_path, stem="native", recovered=False):
    timing = tmp_path / (stem + ".time.txt")
    timing.write_text('Command being timed: "/tool -c 32"\n'
        'User time (seconds): 24.0\nSystem time (seconds): 3.0\n'
        'Elapsed (wall clock) time (h:mm:ss or m:ss): 0:12.00\n'
        'Maximum resident set size (kbytes): 4096\nExit status: 0\n')
    value = {"job_id": "123", "native_argv": ["/tool", "-c", "32"],
             "exit_code": 0, "timing": provenance.record(timing)}
    if recovered:
        value["command"] = value.pop("native_argv")
        value["native_exit_code"] = value.pop("exit_code")
        value["native_timing"] = value.pop("timing")
    return write_json(tmp_path / (stem + ".json"), value)


def test_inventory_is_exact():
    rows = [{"key": k} for k in provenance.KEYS]
    assert list(provenance.inventory(rows)) == list(provenance.KEYS)
    for bad in [rows[:-1], rows + [rows[0]], rows[:-1] + [{"key": "unknown"}]]:
        with pytest.raises(ValueError):
            provenance.inventory(bad)


@pytest.mark.parametrize("value", [None, True, float("nan"), float("inf"), -1, 2, .51])
def test_score_binding_rejects_wrong_values(value):
    with pytest.raises(ValueError):
        provenance.same_score(value, .5)


def test_reader_requires_unchanged_pin(tmp_path):
    pin = write_json(tmp_path / "source.json", {"scope": "old"})
    reader = provenance.Reader()
    assert reader.read(pin) == {"scope": "old"}
    Path(pin["path"]).write_text('{"scope":"new"}')
    with pytest.raises(ValueError):
        reader.read(pin)


@pytest.mark.parametrize("recovered", [False, True])
def test_timed_resource_scope(tmp_path, recovered):
    scope = "recovered downstream mode4 only; excludes BLAST/BPO" if recovered else "full native inference"
    pin = execution(tmp_path, recovered=recovered)
    result = provenance.timed_execution(provenance.Reader(), pin, scope)
    assert result["measurement"]["elapsed_seconds"] == 12
    assert result["full_inference"] is (not recovered)
    assert "not aggregate" in result["memory_scope"]


@pytest.mark.parametrize("damage", ["failed", "command", "scheduler", "ambiguous", "missing"])
def test_native_binding_refusal(tmp_path, damage):
    pin = execution(tmp_path)
    value = json.loads(Path(pin["path"]).read_text())
    if damage in {"ambiguous", "missing"}:
        native = {"checked_records": [dict(pin, path="/one/execution.json"),
                  dict(pin, path="/different/execution.json")]} if damage == "ambiguous" else {}
        with pytest.raises(ValueError):
            provenance.execution_pin(native)
        return
    if damage == "failed":
        value["exit_code"] = 1
    elif damage == "command":
        value["native_argv"] = ["/different"]
    pin = write_json(Path(pin["path"]), value)
    scheduler = {"JobIDRaw": "wrong", "State": "COMPLETED", "ExitCode": "0:0"} if damage == "scheduler" else None
    with pytest.raises(ValueError):
        provenance.timed_execution(provenance.Reader(), pin, "full native inference", scheduler)


@pytest.fixture
def bundle(tmp_path, monkeypatch):
    output = {"path": "/inherited/output", "bytes": 1, "sha256": "a" * 64}
    scores, ob, qfo, three, parity, inputs, of_rows = [], [], [], [], [], [], []
    sonic = None
    for i, key in enumerate(provenance.KEYS):
        scores.append({"key": key, "label": key, "scores": {"OrthoBench": .5, "ThreeKingdoms": .5,
            "QfO_secondary_mean": .5, **{e: .5 for e in ("GO", "EC", "VGNC", "SwissTrees", "TreeFam-A", "FAS")}},
            "orthobench_retained_evidence": {"prediction_provenance": output}, "three_kingdoms": {"groups": output}})
        if i >= 3:
            scores[-1]["orthobench_retained_evidence"] = {}
            scores[-1]["orthobench_supplemental_readback"] = {"prediction": output}
        ob.append(dict(key=key, prediction=output, weighted_refog_f1=.5, input_note="Unknown",
            inference_wall_seconds=None, conversion_wall_seconds=None, peak_process_rss_kib=None, timing_basis="Unknown", output_semantics="groups"))
        counts = {"f_score": .5}
        three.append(dict(key=key, groups=output, counts=counts, input_status="Unknown", semantics="BUSCO co-membership", run="historical", use="comparison"))
        parity.append(dict(key=key, version="3.1.5" if key.startswith("orthofinder") else "declared",
            provenance={"orthogroups_sha256": output["sha256"]}, score=counts,
            performance={"runtime_kind": "measured", "wall_s": 99, "peak_rss_kib": 77,
                         "memory_measurement": "old scope", "cpus_requested": 32}))
        inputs.append(dict(method=key, evidence=[], copies=[]))
        recovered = key == "orthomcl_1_4"
        pin = execution(tmp_path, f"exec{i}", recovered)
        native = {"execution": pin, "native_pairs": output}
        native_pin = write_json(tmp_path / f"admission{i}.json", native)
        conversion = dict(semantics="pairs", participant="participant", total_pairs=2, retained_pairs=2,
            removed_mapping_pairs=0, filtered_pairs=output, admission=native_pin)
        if key == "orthohmm_phylogeny_satellite_v2":
            conversion.pop("admission")
            conversion.update(native_admission=native_pin, native_input=output)
        conversion_pin = write_json(tmp_path / f"conversion{i}.json", conversion)
        endpoints = {e: .5 for e in ("GO", "EC", "VGNC", "SwissTrees", "TreeFam-A", "FAS")}
        qfo.append(dict(key=key, conversion=conversion_pin, prediction_semantics="pairs", submitted_pairs=2,
            retained_pairs=2, removed_mapping_pairs=0, status="admitted", scores=endpoints, secondary_mean=.5,
            admission=native_pin, details={}))
        if key.startswith("orthofinder"):
            of_rows.append(dict(key=key, pair_file=output, native_admission=native_pin, scores=endpoints, conversion_wall_seconds=1))
        if key == "sonicparanoid_2_0_9":
            sonic = dict(normalized=output, counts=counts, execution=pin,
                scheduler={"JobIDRaw": "123", "State": "COMPLETED", "ExitCode": "0:0"})
    docs = {"scores": {"rows": scores}, "ob": {"rows": ob}, "qfo": {"methods": qfo}, "three": {"rows": three},
        "parity": {"methods": parity}, "sonic": sonic, "three_inputs": {"methods": inputs}, "mcl_inputs": {},
        "of": {"rows": of_rows, "native_command": ["/OF"], "full_native_resources": {"elapsed_seconds": 100},
               "native_scheduler": {}, "checked_records": []}}
    base = tmp_path / "repo" / "benchmark_tools" / "results"
    base.mkdir(parents=True)
    sources = {}
    for name, doc in docs.items():
        pin = write_json(base / (name + ".json"), doc)
        sources[name] = (name + ".json", pin["sha256"])
    monkeypatch.setattr(provenance, "SOURCES", sources)
    return tmp_path / "repo", docs, sources


def test_complete_24_row_consolidation_preserves_scopes(bundle, tmp_path):
    repo, _, _ = bundle
    result = provenance.assemble(repo, tmp_path / "register")
    assert len(result["rows"]) == 24
    assert len({(r["dataset"], r["key"]) for r in result["rows"]}) == 24
    assert result["complete_transitive_provenance"] is False
    assert result["publication_ready"] is False
    sonic = next(r for r in result["rows"] if r["dataset"] == "ThreeKingdoms" and r["key"].startswith("sonic"))
    assert sonic["resources"][0]["measurement"]["elapsed_seconds"] == 12  # Never historical99.
    of_seq = next(r for r in result["rows"] if r["dataset"] == "QfO" and r["key"].endswith("sequence_only"))
    assert all(r["full_inference"] is False for r in of_seq["resources"])
    fastoma = next(r for r in result["rows"] if r["dataset"] == "QfO" and r["key"].startswith("fastoma"))
    assert "Docker" in fastoma["resources"][0]["memory_scope"]
    assert (tmp_path / "register" / "resources.tsv").exists()


@pytest.mark.parametrize("damage", ["ob_score", "ob_output", "three_output", "sonic_old_output", "qfo_score", "native_binding", "changed_source"])
def test_consolidation_refuses_mixed_runs(bundle, tmp_path, monkeypatch, damage):
    repo, docs, sources = bundle
    if damage == "ob_score":
        docs["ob"]["rows"][0]["weighted_refog_f1"] = .4
    elif damage == "ob_output":
        docs["ob"]["rows"][0]["prediction"] = {"sha256": "wrong"}
    elif damage == "three_output":
        docs["parity"]["methods"][0]["provenance"]["orthogroups_sha256"] = "wrong"
    elif damage == "sonic_old_output":
        docs["sonic"]["normalized"] = {"sha256": "historical-output"}
    elif damage == "qfo_score":
        docs["qfo"]["methods"][0]["scores"]["GO"] = .4
    elif damage == "native_binding":
        docs["of"]["rows"][0]["pair_file"] = {"sha256": "different"}
    sources = copy.deepcopy(sources)
    for name, doc in docs.items():
        pin = write_json(repo / "benchmark_tools/results" / (name + ".json"), doc)
        if damage != "changed_source":
            sources[name] = (name + ".json", pin["sha256"])
    if damage == "changed_source":
        (repo / "benchmark_tools/results/scores.json").write_text("{}")
    monkeypatch.setattr(provenance, "SOURCES", sources)
    with pytest.raises((ValueError, KeyError)):
        provenance.assemble(repo, tmp_path / "rejected")
    assert not (tmp_path / "rejected").exists()


def test_existing_destination_is_not_overwritten(bundle, tmp_path):
    repo, _, _ = bundle
    output = tmp_path / "existing"
    output.mkdir()
    marker = output / "marker"
    marker.write_text("preserve")
    with pytest.raises(FileExistsError):
        provenance.assemble(repo, output)
    assert marker.read_text() == "preserve"


@pytest.mark.parametrize("damage", ["missing", "conflicting"])
def test_selected_ob_prediction_requires_bound_evidence(damage):
    row = {"orthobench_retained_evidence": {}}
    if damage == "conflicting":
        row["orthobench_retained_evidence"]["prediction_provenance"] = {"sha256": "first"}
        row["orthobench_supplemental_readback"] = {"prediction": {"sha256": "second"}}
    with pytest.raises(ValueError):
        provenance.ob_prediction(row)
