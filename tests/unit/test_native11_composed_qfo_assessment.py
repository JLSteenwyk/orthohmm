"""Real endpoint/trace/FAS admission kernels; scheduler/native contexts are stubs."""

from copy import deepcopy
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import pytest

from benchmark_tools import run_native11_composed_qfo_assessment as runner
from benchmark_tools import admit_native11_composed_qfo_assessment as admission
from benchmark_tools import run_allocated_native_factorial_qfo_assessment as ordinary
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_run_native_factorial_qfo_assessment import (
    joined as original_joined, short_root, put, write_native_endpoints, stage_fixture)


def stage():
    value, producer = stage_fixture(index=11)
    value.update(schema=runner.converter.SCHEMA, status=runner.converter.STATUS, native_job_id=23985,
        original_review_translated=False, conversion_started_monotonic_ns=1, conversion_finished_monotonic_ns=2,
        composed_binding=dict(composed_schema_preserved=True, original_review_translated=False, next_identity_authorized=False))
    return value, producer


def test_explicit_stage_is_not_ordinary_stage():
    value, producer = stage()
    runner.validate_stage(value, producer)
    with pytest.raises(ValueError):
        ordinary.validate_stage(value, producer)


@pytest.mark.parametrize("field", ["schema", "status", "native_index", "native_job_id", "cell", "participant",
    "conversion_kind", "semantics", "accuracy_evaluated", "original_review_translated", "next_identity_authorized",
    "automatic_retry", "total_pairs", "retained_pairs", "expected_pairs", "removed_mapping_pairs", "empty_predictions",
    "composed_binding", "pair_coverage", "conversion_started_monotonic_ns"])
def test_changed_stage_identity_semantics_counts_or_proof_is_rejected(field):
    value, producer = stage()
    value[field] = not value[field] if type(value[field]) is bool else "wrong"
    with pytest.raises((ValueError, TypeError, AttributeError)):
        runner.validate_stage(value, producer)


@pytest.fixture
def joined(original_joined, monkeypatch):
    root, original_ref, original_stage, manifest = original_joined
    tools = root / "benchmark_tools"
    original = Path(runner.__file__).parent
    names = ("native11_composed_review_binding.py", "prepare_native11_composed_qfo_pairs.py",
        "run_native11_composed_qfo_assessment.py", "admit_native11_composed_qfo_assessment.py",
        "prepare_native_factorial_qfo_pairs.py", "admit_qfo_recovered_assessment.py", "audit_qfo_fas_samples.py",
        "validate_qfo_native_assessment.py", "qfo_summarize_scores.py", "admit_native_factorial_qfo_assessment.py")
    for name in names:
        shutil.copyfile(original / name, tools / name)
    for module in (runner, runner.converter, admission):
        monkeypatch.setattr(module, "__file__", str(tools / Path(module.__file__).name))
        monkeypatch.setattr(module, "ROOT", root)
    monkeypatch.setattr(runner.converter, "DESTINATION", root / "composed_conversion")
    monkeypatch.setattr(admission, "DESTINATION", root / "composed_admission")
    value, producer = stage()
    value.update({k: deepcopy(original_stage[k]) for k in ("pairs", "filtered_pairs", "input_fastas", "mapping",
        "environment_manifest", "gene_ownership_sha256", "request", "plan", "terminal_review")})
    value["source"] = record(runner.converter.__file__)
    value["binding_source"] = record(tools / "native11_composed_review_binding.py")
    value["conversion_kernel_source"] = record(tools / "prepare_native_factorial_qfo_pairs.py")
    value["conversion_log"] = put(root / "composed_conversion/conversion.log", "fixture group conversion\n")
    output_root = root / "attempt"
    value["native_input"] = put(output_root / "native/orthohmm_working_res/orthohmm_edges_clustered.txt", "fixture groups\n")
    value["allocated_ready"] = put(output_root / "measurement/ready.json", dict(fixture=True))
    value["native_cpu_ids"] = [1]
    value["amendment"] = put(root / "amendment.json", {})
    value["checked_records"] = [value["source"], value["native_input"]]
    ref = put(root / "composed_conversion/results.json", value)
    request = dict(plan=value["plan"], amendment=value["amendment"], job_id=23985)
    run = dict(index=11, cell=value["cell"], genes=3, inputs=value["input_fastas"], output_root=str(output_root))
    outputs = dict(native_cpu_ids=[1], allocated_ready=value["allocated_ready"],
        gene_ownership_sha256=value["gene_ownership_sha256"], checked_files=[value["native_input"]])
    context = [request, {}, {}, run, {}, outputs, "group", {"fixture_native_scheduler": True}, [],
               deepcopy(value["composed_binding"])]
    monkeypatch.setattr(runner, "native_binding", lambda *a: tuple(context))
    monkeypatch.setattr(runner.converter, "ENV_SHA", value["environment_manifest"]["sha256"])
    monkeypatch.setattr(runner, "accounting", lambda job: ("fixture conversion", deepcopy(producer)))
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setenv("SLURM_JOB_ID", "778")
    return root, ref, value, manifest


def execute(joined, monkeypatch, *, exit_code=0, tamper=False):
    root, ref, value, manifest = joined
    def run(command, **kwargs):
        assert command[command.index("--challenges_ids") + 1] == "GO EC VGNC SwissTrees TreeFam-A FAS"
        assert command[command.index("--event_year") + 1] == "2020"
        assert "-resume" not in command
        assert kwargs["env"]["NXF_OFFLINE"] == "true"
        write_native_endpoints(Path(command[command.index("--results_dir") + 1]), value, manifest)
        if tamper:
            Path(value["filtered_pairs"]["path"]).write_text("changed")
        return SimpleNamespace(returncode=exit_code)
    monkeypatch.setattr(runner.subprocess, "run", run)
    return runner.run(root, ref, "777", record(runner.__file__)["sha256"])


def admit_fixture(joined, monkeypatch):
    root, ref, value, manifest = joined
    monkeypatch.setattr(admission, "accounting", lambda job: ("fixture assessment", dict(
        JobIDRaw="778", State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="8")))
    return admission.admit(root, ref, "777", "778", admission.DESTINATION, record(admission.__file__)["sha256"])


def test_read_only_spec_uses_same_endpoints_and_new_paths(joined):
    root, ref, value, manifest = joined
    spec = runner.run(root, ref, "777", record(runner.__file__)["sha256"], check_only=True)
    assert spec["schema"] == runner.SCHEMA and spec["status"] == "prepared_unrun"
    assert spec["work"].endswith("/cn11")
    assert spec["results"].endswith("/composed_native11_v1")
    assert not Path(spec["cwd"]).exists()
    assert spec["accuracy_admitted"] is False


def test_real_independent_endpoint_trace_and_fas_admission(joined, monkeypatch):
    execution = execute(joined, monkeypatch)
    assert execution["status"] == "process_succeeded_pending_independent_admission"
    assert execution["accuracy_admitted"] is False
    result = admit_fixture(joined, monkeypatch)
    assert result["schema"] == admission.SCHEMA
    assert result["status"] == "native11_composed_qfo_assessment_admitted"
    assert result["accuracy_admitted"] is True
    assert result["publication_ready"] is result["next_identity_authorized"] is False
    assert result["original_review_translated"] is False
    assert result["assessment"]["participant"] == runner.converter.PARTICIPANT


@pytest.mark.parametrize("problem", ["failure", "pair_tamper"])
def test_failed_or_tampered_scoring_retains_failure_without_admission(joined, monkeypatch, problem):
    with pytest.raises(ValueError):
        execute(joined, monkeypatch, exit_code=1 if problem == "failure" else 0, tamper=problem == "pair_tamper")
    report = json.loads((joined[0] / "benchmarks/results/native11_composed_qfo_assessment_v1/results.json").read_text())
    assert report["status"] == "failed"
    assert report["accuracy_admitted"] is False


@pytest.mark.parametrize("problem", ["endpoint", "trace", "fas", "inventory"])
def test_invalid_real_endpoint_or_inventory_cannot_be_admitted(joined, monkeypatch, problem):
    execution = execute(joined, monkeypatch)
    results = Path(execution["results"])
    if problem == "endpoint":
        path = results / "assessment_out/Assessment_datasets.json"
        values = json.loads(path.read_text())
        values[0]["participant_id"] = "wrong"
        put(path, values)
    elif problem == "trace":
        path = results / "stats/trace_fixture.txt"
        path.write_text(path.read_text().replace("COMPLETED", "FAILED", 1))
    elif problem == "fas":
        import gzip
        with gzip.open(results / "results/FAS/sample_raw.txt.gz", "wt") as stream:
            stream.write("Acc1\tAcc2\tFAS\nA1\tB1\t0.9\nA1\tC1\t0.9\n")
    else:
        put(results / "extra.txt", "unrecorded")
    # Refreshing the inventory deliberately exercises semantic rejection beyond SHA checks.
    if problem != "inventory":
        execution["outputs"] = [record(p) for p in sorted(results.rglob("*")) if p.is_file()]
        put(Path(execution["cwd"]) / "results.json", execution)
    with pytest.raises(ValueError):
        admit_fixture(joined, monkeypatch)
    report = json.loads((admission.DESTINATION / "results.json").read_text())
    assert report["accuracy_admitted"] is False
    assert report["status"] == "native11_composed_qfo_admission_failed_retained"
