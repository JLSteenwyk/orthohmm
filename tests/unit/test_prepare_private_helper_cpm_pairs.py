import copy
import gzip
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import prepare_private_helper_cpm_pairs as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


JOB = "24001"


def accounting(state="COMPLETED", exit_code="0:0", node="bizon", cpus="2", mem="64G"):
    header = "JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS|ReqMem\n"
    return (header + f"{JOB}|{JOB}|{state}|{exit_code}|00:01:00|{node}|{cpus}|{mem}\n"
            + f"{JOB}.batch|{JOB}.batch|COMPLETED|0:0|00:01:00|bizon|2|\n")


@pytest.mark.parametrize("change", [{"state": "RUNNING"}, {"state": "PENDING"},
    {"state": "FAILED"}, {"state": "TIMEOUT"}, {"state": "CANCELLED"},
    {"exit_code": "1:0"}, {"node": "other"}, {"cpus": "32"}, {"mem": "192G"}])
def test_terminal_admission_required(change):
    with pytest.raises(ValueError):
        module.completed_admission(accounting(**change), JOB)


@pytest.mark.parametrize("job", [None, True, 24001, "", "0", "024001", "24001_1", "24001.batch", "-1", "１２"])
def test_standalone_job_identity(job):
    with pytest.raises(ValueError):
        module.completed_admission("", job)


def test_parent_and_steps_unique_complete():
    text = accounting()
    assert module.completed_admission(text, JOB)["AllocCPUS"] == "2"
    lines = text.splitlines()
    for broken in ("\n".join(lines[:2]), text + lines[1] + "\n", text + lines[2] + "\n",
                   text.replace("24001.batch|24001.batch|COMPLETED|0:0", "24001.batch|24001.batch|FAILED|1:0"),
                   text.replace("24001|24001|", "24001|24002|")):
        with pytest.raises(ValueError):
            module.completed_admission(broken, JOB)


def scheduled(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(module.sys, "executable", module.PYTHON)
    for key, value in {"SLURM_JOB_ID": "24002", "SLURM_CPUS_PER_TASK": "2",
        "SLURM_JOB_NODELIST": "bizon", "SLURM_MEM_PER_NODE": "65536"}.items():
        monkeypatch.setenv(key, value)
    monkeypatch.delenv("SLURM_ARRAY_TASK_ID", raising=False)


@pytest.mark.parametrize("key,value", [("SLURM_JOB_ID", ""), ("SLURM_CPUS_PER_TASK", "32"),
    ("SLURM_JOB_NODELIST", "other"), ("SLURM_MEM_PER_NODE", "196608"), ("SLURM_ARRAY_TASK_ID", "1")])
def test_allocation_before_input_reads(tmp_path, monkeypatch, key, value):
    scheduled(tmp_path, monkeypatch)
    monkeypatch.setenv(key, value)
    monkeypatch.setattr(module, "accounting", lambda *_: pytest.fail("Read scheduler before allocation gate"))
    with pytest.raises(ValueError, match="scheduled"):
        module.prepare(tmp_path, JOB, "", "", "")


def test_wrong_controller_or_cwd(tmp_path, monkeypatch):
    scheduled(tmp_path, monkeypatch)
    monkeypatch.setattr(module.sys, "executable", "/unreviewed/python")
    with pytest.raises(ValueError, match="controller"):
        module.prepare(tmp_path, JOB, "", "", "")
    monkeypatch.setattr(module.sys, "executable", module.PYTHON)
    with pytest.raises(ValueError, match="verification directory"):
        module.prepare(tmp_path / "other", JOB, "", "", "")


def setup(tmp_path, monkeypatch):
    scheduled(tmp_path, monkeypatch)
    executor = tmp_path / module.ADMITTER_EXECUTOR
    filenames = ("benchmark_tools/admit_private_helper_cpm_phylogeny.py", module.native_checker.PROTOCOL,
        "benchmark_tools/results/qfo_private_cpm_native_admission_20261001.sh",
        "tests/unit/test_admit_private_helper_cpm_phylogeny.py")
    root = Path(module.__file__).resolve().parents[1]
    for name in filenames:
        target = executor / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((root / name).read_bytes())
    protocol = tmp_path / module.native_checker.PROTOCOL
    protocol.parent.mkdir(parents=True, exist_ok=True)
    protocol.write_bytes((root / module.native_checker.PROTOCOL).read_bytes())
    parent = tmp_path / module.native_checker.SUBMISSION
    parent.write_text("{}")
    monkeypatch.setattr(module, "NATIVE_SUBMISSION_SHA", record(parent)["sha256"])
    submitted = dict(status="private_recovered_qfo_high_cpm_native_admission_submitted", job_id=JOB,
        executor=str(executor), executor_commit=module.ADMITTER_COMMIT, executor_clean=True,
        native_inference=False, accuracy_evaluated=False, controlled_timing=False, publication_ready=False,
        submission_argv=["sbatch", "--parsable", str(executor / filenames[2]), str(executor),
                         module.ADMITTER_COMMIT, module.NATIVE_PROTOCOL_SHA],
        parent_submission=record(parent), source_records=[record(executor / name) for name in filenames])
    submission = tmp_path / f"benchmark_tools/results/qfo_private_cpm_native_admission_submission_{JOB}.json"
    submission.write_text(json.dumps(submitted))
    converter_protocol = tmp_path / module.PROTOCOL
    converter_protocol.write_text("# Prospective conversion fixture\n")
    pairs = tmp_path / "native.tsv"
    pairs.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\nA\ts1\tB\ts2\nA\ts1\tC\ts3\n")
    metadata = tmp_path / "metadata.json"
    metadata.write_text("{}")
    candidates = tmp_path / "recovered_candidates.txt"
    candidates.write_text("A B C\n")
    candidate_admission = tmp_path / "candidate_admission.json"
    private_admission = tmp_path / "private_admission.json"
    candidate_admission.write_text("{}")
    private_admission.write_text("{}")
    verified = dict(baseline=dict(manifest=dict(input_fastas=[])),
        candidates=dict(admission_record=record(candidate_admission), arm=dict(partition=record(candidates))),
        private_control=dict(admission=record(private_admission)), checked_records=[record(candidates)])
    cell = dict(label="candidate_cpm_high", argv=["private-python", "inferred-phylogeny"])
    native = dict(status="private_recovered_cpm_native_pairs_verified_unscored", arm="cpm_high", index=1,
        source=submitted["source_records"][0], cell=cell, native_pair_count=2, native_pairs=record(pairs),
        candidate_admission=verified["candidates"]["admission_record"],
        private_native_admission=verified["private_control"]["admission"],
        protocol=record(protocol), submission=record(parent),
        native_group_integrity=dict(native_manifest=record(metadata)), accuracy_evaluated=False,
        scoring_admitted=False, controlled_timing=False, publication_ready=False)
    native["checked_records"] = [native[key] for key in ("source", "protocol", "submission", "native_pairs")]
    native["checked_records"].append(record(metadata))
    admission = tmp_path / module.ADMISSION
    admission.parent.mkdir(parents=True, exist_ok=True)
    admission.write_text(json.dumps(native))
    mapping = tmp_path / "mapping.json.gz"
    with gzip.open(mapping, "wt") as stream:
        json.dump({"mapping": dict.fromkeys("ABC", 1)}, stream)
    environment = tmp_path / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment.write_text(json.dumps(dict(reference_files=[record(mapping)])))
    monkeypatch.setattr(module, "ENV_SHA", record(environment)["sha256"])
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: module.ADMITTER_COMMIT)
    monkeypatch.setattr(module, "accounting", lambda *_: accounting())
    checks = []
    def verify(root):
        checks.append(root)
        return verified
    producer = SimpleNamespace(OUTPUT="native-output", verify_sources=verify,
        native_command=lambda *args: (cell, None, None, None))
    monkeypatch.setattr(module.native_checker, "producer", lambda *args: (producer, verify(tmp_path)))
    monkeypatch.setattr(module, "gene_ownership", lambda *args: ({"A": "s1", "B": "s2", "C": "s3"}, {}))
    def run(command, **kwargs):
        if command[0] == "git":
            return
        assert command[:3] == [module.PYTHON, "-B", native["source"]["path"]]
        assert command[command.index("--submission-sha256") + 1] == module.NATIVE_SUBMISSION_SHA
        assert command[command.index("--protocol-sha256") + 1] == module.NATIVE_PROTOCOL_SHA
        Path(command[command.index("--output") + 1]).write_text(json.dumps(native))
    monkeypatch.setattr(module.subprocess, "run", run)
    return SimpleNamespace(native=native, verified=verified, cell=cell, submitted=submitted,
        admission=admission, submission=submission, protocol=converter_protocol, mapping=mapping,
        output=tmp_path / module.OUTPUT, checks=checks, producer=producer, pairs=pairs, metadata=metadata)


def invoke(root, data):
    return module.prepare(root, JOB, record(data.admission)["sha256"],
                          record(data.submission)["sha256"], record(data.protocol)["sha256"])


def test_full_lossless_prepare_unchanged_helpers(tmp_path, monkeypatch):
    data = setup(tmp_path, monkeypatch)
    value = invoke(tmp_path, data)
    assert value["status"] == "private_recovered_cpm_native_pairs_prepared_unscored"
    assert value["arm"] == "cpm_high" and value["index"] == 1
    assert value["participant"] == "ohmm_qfo_parameter_cpm_high"
    assert value["semantics"] == "native phylogenetically inferred pairs"
    assert value["written_pairs"] == value["total_pairs"] == value["retained_pairs"] == 2
    assert value["removed_mapping_pairs"] == 0
    assert all(value[key] is False for key in
               ("accuracy_evaluated", "scoring_admitted", "controlled_timing", "publication_ready"))
    assert (data.output / "pairs.qfo.tsv").read_text() == "A\tB\nA\tC\n"
    assert value["pairs"]["sha256"] == value["filtered_pairs"]["sha256"]
    assert value["native_admission_recheck"] in value["checked_records"]
    assert len(data.checks) == 2
    assert not (data.output / "pairs.partial.tsv").exists()
    with pytest.raises(FileExistsError):
        invoke(tmp_path, data)


@pytest.mark.parametrize("key,value", [("status", "old_cpm_admission"), ("arm", "cpm_low"),
    ("index", True), ("source", {}), ("cell", {}), ("candidate_admission", {}),
    ("private_native_admission", {}), ("accuracy_evaluated", True), ("scoring_admitted", True),
    ("controlled_timing", 0), ("publication_ready", 0), ("native_pair_count", True),
    ("native_pair_count", 0), ("checked_records", [])])
def test_wrong_native_arm_provenance_or_claim_rejected(tmp_path, monkeypatch, key, value):
    data = setup(tmp_path, monkeypatch)
    native = copy.deepcopy(data.native)
    native[key] = value
    with pytest.raises(ValueError):
        module.validate_native(native, data.native["source"], data.verified, data.cell, tmp_path)


@pytest.mark.parametrize("key,value", [("status", "wrong"), ("job_id", "24003"),
    ("executor_commit", "wrong"), ("executor_clean", False), ("native_inference", True),
    ("accuracy_evaluated", 0), ("source_records", []), ("submission_argv", []), ("parent_submission", {})])
def test_changed_submission_rejected_before_outputs(tmp_path, monkeypatch, key, value):
    data = setup(tmp_path, monkeypatch)
    data.submitted[key] = value
    data.submission.write_text(json.dumps(data.submitted))
    with pytest.raises(ValueError):
        invoke(tmp_path, data)
    assert not data.output.exists()


def test_active_admission_precedes_data_reads(tmp_path, monkeypatch):
    scheduled(tmp_path, monkeypatch)
    monkeypatch.setattr(module, "accounting", lambda *_: accounting(state="RUNNING"))
    monkeypatch.setattr(module, "record", lambda *_: pytest.fail("Read data while admission live"))
    with pytest.raises(ValueError, match="completed"):
        module.prepare(tmp_path, JOB, "", "", "")


@pytest.mark.parametrize("problem", ["fresh", "fresh_failure", "mapping", "runtime", "count",
    "pair_order", "species", "metadata", "partial_bytes", "final_accounting", "final_evidence"])
def test_failed_conversion_preserves_unscored_evidence(tmp_path, monkeypatch, problem):
    data = setup(tmp_path, monkeypatch)
    if problem == "fresh":
        def run(command, **kwargs):
            if command[0] != "git":
                Path(command[command.index("--output") + 1]).write_text("{}")
        monkeypatch.setattr(module.subprocess, "run", run)
    elif problem == "fresh_failure":
        def run(command, **kwargs):
            if command[0] != "git":
                raise module.subprocess.CalledProcessError(1, command)
        monkeypatch.setattr(module.subprocess, "run", run)
    elif problem == "mapping":
        with gzip.open(data.mapping, "wt") as stream:
            json.dump({"mapping": dict.fromkeys("AB", 1)}, stream)
        environment = tmp_path / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
        environment.write_text(json.dumps(dict(reference_files=[record(data.mapping)])))
        monkeypatch.setattr(module, "ENV_SHA", record(environment)["sha256"])
    elif problem == "runtime":
        data.producer.verify_sources = lambda *_: {}
    elif problem in {"count", "pair_order", "species"}:
        if problem == "count":
            data.native["native_pair_count"] = 3
        else:
            text = data.pairs.read_text()
            if problem == "pair_order":
                text = text.replace("A\ts1\tC\ts3", "A\ts1\tB\ts2")
            else:
                text = text.replace("B\ts2", "B\ts1")
            data.pairs.write_text(text)
            old = data.native["native_pairs"]
            data.native["native_pairs"] = record(data.pairs)
            data.native["checked_records"].remove(old)
            data.native["checked_records"].append(data.native["native_pairs"])
        data.admission.write_text(json.dumps(data.native))
    elif problem == "metadata":
        monkeypatch.setattr(module, "gene_ownership", lambda *_: (_ for _ in ()).throw(ValueError("Bad universe")))
    elif problem == "partial_bytes":
        original = module.convert
        def convert(*args):
            value = original(*args)
            (data.output / "pairs.qfo.partial.tsv").write_text("A\tB\n")
            return value
        monkeypatch.setattr(module, "convert", convert)
    elif problem == "final_accounting":
        observed = []
        def later(*args):
            observed.append(args)
            return accounting(state="FAILED" if len(observed) > 1 else "COMPLETED")
        monkeypatch.setattr(module, "accounting", later)
    else:
        def verify(*args):
            data.metadata.write_text("{\"changed\":true}")
            return data.verified
        data.producer.verify_sources = verify
    with pytest.raises((ValueError, module.subprocess.CalledProcessError)):
        invoke(tmp_path, data)
    assert (data.output / "preflight.json").exists()
    saved = json.loads((data.output / "results.json").read_text())
    assert saved["status"] == "failed"
    assert all(saved[key] is False for key in
               ("accuracy_evaluated", "scoring_admitted", "controlled_timing", "publication_ready"))
    assert not (data.output / "pairs.tsv").exists()
    if problem == "mapping":
        assert json.loads((data.output / "conversion_counts.json").read_text())["removed_mapping_pairs"] == 1


def test_changed_admission_pin_rejected_before_outputs(tmp_path, monkeypatch):
    data = setup(tmp_path, monkeypatch)
    with pytest.raises(ValueError):
        module.prepare(tmp_path, JOB, "0" * 64, record(data.submission)["sha256"], record(data.protocol)["sha256"])
    assert not data.output.exists()


def test_symlink_output_rejected_before_scheduler(tmp_path, monkeypatch):
    scheduled(tmp_path, monkeypatch)
    output = tmp_path / module.OUTPUT
    output.parent.mkdir(parents=True)
    output.symlink_to(tmp_path / "missing", target_is_directory=True)
    monkeypatch.setattr(module, "accounting", lambda *_: pytest.fail("Read scheduler after output collision"))
    with pytest.raises(FileExistsError):
        module.prepare(tmp_path, JOB, "", "", "")


@pytest.mark.parametrize("problem", ["revision", "source", "dirty"])
def test_source_revision_or_dirty_executor_rejected(tmp_path, monkeypatch, problem):
    data = setup(tmp_path, monkeypatch)
    if problem == "revision":
        monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "changed-revision")
    elif problem == "source":
        checker = tmp_path / module.ADMITTER_EXECUTOR / "benchmark_tools/admit_private_helper_cpm_phylogeny.py"
        checker.write_text("# changed\n")
        data.submitted["source_records"][0] = record(checker)
        data.submission.write_text(json.dumps(data.submitted))
    else:
        def run(command, **kwargs):
            raise module.subprocess.CalledProcessError(1, command)
        monkeypatch.setattr(module.subprocess, "run", run)
    with pytest.raises((ValueError, module.subprocess.CalledProcessError)):
        invoke(tmp_path, data)
    assert not data.output.exists()


def test_ambiguous_reference_mapping_rejected(tmp_path, monkeypatch):
    data = setup(tmp_path, monkeypatch)
    environment = tmp_path / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment.write_text(json.dumps(dict(reference_files=[record(data.mapping)] * 2)))
    monkeypatch.setattr(module, "ENV_SHA", record(environment)["sha256"])
    with pytest.raises(ValueError, match="one frozen"):
        invoke(tmp_path, data)
    assert not data.output.exists()


def test_normalization_collision_preserved_as_failure(tmp_path, monkeypatch):
    data = setup(tmp_path, monkeypatch)
    monkeypatch.setattr(module, "gene_ownership", lambda *a: ({"sp|A|one": "s1", "A": "s2"}, {}))
    with pytest.raises(ValueError, match="normalized accession"):
        invoke(tmp_path, data)
    assert json.loads((data.output / "results.json").read_text())["status"] == "failed"
    assert not (data.output / "pairs.tsv").exists()


def test_launch_script_has_fixed_resources_and_explicit_bindings():
    script = Path(module.__file__).parent / "results/qfo_private_cpm_pair_conversion_20261001.sh"
    text = script.read_text()
    for required in ("--nodelist=bizon", "--cpus-per-task=2", "--mem=64G", "--time=04:00:00",
        "--no-requeue", "prepare_private_helper_cpm_pairs.py", "--admission-job", "--admission-sha256",
        "--admission-submission-sha256", "--protocol-sha256", "/home/bizon/anaconda3/bin/python -B"):
        assert required in text
    for forbidden in ("--array", "--dependency", "--gres", "--exclusive", "sbatch ", "while ", "until "):
        assert forbidden not in text
