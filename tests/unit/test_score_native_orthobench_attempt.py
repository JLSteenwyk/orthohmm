"""Terminal scoring gates and frozen RefOG semantics; no native inference."""

import itertools
import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools import score_native_orthobench_attempt as score
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.validate_native_factorial_outputs import Evidence


def admission_fixture():
    request_ref = dict(path="/request.json", bytes=1, sha256="a" * 64)
    plan_ref = dict(path="/plan.json", bytes=1, sha256="b" * 64)
    request = dict(plan=plan_ref, job_id=22428)
    run = dict(dataset="orthobench", index=1, cell="p0_c0_r1", repeat=0)
    review = dict(schema="native_factorial_terminal_review_v1", request=request_ref, plan=plan_ref,
        job_id=22428, index=1, cell="p0_c0_r1", dataset="orthobench", repeat=0, status="native_success",
        scheduler_state="COMPLETED", scheduler_exit_code="0:0", terminal_reviewed=True,
        primary_resources_replayed=True, shared_host_resources_reviewed=True, native_outputs_validated=True,
        execution_scope=score.SCOPE, resource_scopes=score.SCOPES, automatic_retry=False, uncontended_timing=False,
        reviews=dict.fromkeys(("runtime", "resources", "environment", "outputs_or_failure")))
    return review, request_ref, request, run


def test_bound_successful_terminal_admission():
    score.admit_score(*admission_fixture())


@pytest.mark.parametrize("key,value", [("job_id", 999), ("index", 0), ("cell", "p0_c1_r1"),
    ("dataset", "qfo_corrected"), ("repeat", 1), ("status", "native_failure_retained"),
    ("schema", "live_start"), ("scheduler_state", "RUNNING"), ("scheduler_exit_code", "1:0"),
    ("terminal_reviewed", False), ("primary_resources_replayed", False),
    ("shared_host_resources_reviewed", False), ("native_outputs_validated", False),
    ("execution_scope", "isolated"), ("resource_scopes", {}), ("automatic_retry", True),
    ("uncontended_timing", True), ("request", {}), ("plan", {}), ("reviews", {})])
def test_pending_failed_tampered_reviews_not_scores(key, value):
    review, ref, request, run = admission_fixture()
    review[key] = value
    with pytest.raises(ValueError, match="successful terminal"):
        score.admit_score(review, ref, request, run)


def test_qfo_is_not_a_refog_endpoint():
    review, ref, request, run = admission_fixture()
    run["dataset"] = "qfo_corrected"
    with pytest.raises(ValueError, match="separate native endpoint"):
        score.admit_score(review, ref, request, run)


def partitions_fixture():
    originals = {frozenset("ab"), frozenset("cdef"), frozenset("u")}
    current = {frozenset("abu"), frozenset("cd"), frozenset("e"), frozenset("f")}
    references = dict(RefOG001=set("ab"), RefOG002=set("cdef"))
    uncertain = dict(RefOG001=set("u"))
    old = (set("abcdefu"), originals)
    new = (set("abcdefu"), current)
    old_point = score.score_partition(originals, references, uncertain)
    return new, old, references, uncertain, old_point


def test_weighted_statistic_not_mean_of_family_f1_and_low_certainty_excluded():
    result = score.scored_partitions(*partitions_fixture())
    assert result["score_fraction"]["f_score"] == pytest.approx(8 / 13)
    assert result["score_fraction"]["f_score"] != pytest.approx((1 + 2 / 7) / 2)
    assert result["score_fraction"]["precision"] == 1
    assert result["difference_from_original_percentage_points"]["f_score"] < 0
    assert result["canonical_partition_comparison"]["partition_equal"] is False
    assert "native_only_groups" in result["canonical_partition_comparison"]
    assert "scaling_only_groups" not in result["canonical_partition_comparison"]
    assert result["reference_genes"] == result["covered_reference_genes"] == 6


def test_identical_partitions_are_label_independent():
    new, old, references, uncertain, old_point = partitions_fixture()
    result = score.scored_partitions(old, old, references, uncertain, old_point)
    assert result["canonical_partition_comparison"]["partition_equal"] is True
    assert all(value == 0 for value in result["difference_from_original_percentage_points"].values())


@pytest.mark.parametrize("value", [float("nan"), float("inf"), True, -1, 50])
def test_cached_point_does_not_override_actual_scoring(value):
    new, old, references, uncertain, old_point = partitions_fixture()
    old_point["f_score"] = value
    with pytest.raises(ValueError, match="does not reproduce"):
        score.scored_partitions(new, old, references, uncertain, old_point)


def test_universe_mismatch_refused():
    new, old, references, uncertain, old_point = partitions_fixture()
    with pytest.raises(ValueError, match="universe differs"):
        score.scored_partitions((new[0] - {"u"}, new[1]), old, references, uncertain, old_point)
    references["RefOG003"] = set("xy")
    with pytest.raises(ValueError, match="RefOG genes missing"):
        score.scored_partitions(new, old, references, uncertain, old_point)


@pytest.fixture
def frozen_fixture(tmp_path, monkeypatch):
    refs = tmp_path / "references"
    refs.mkdir()
    uncertain = refs / "low_certainty_assignments"
    uncertain.mkdir()
    records, references = [], {}
    lines = []
    for i in range(1, 71):
        genes = [f"g{i:03d}a", f"g{i:03d}b"]
        name = f"RefOG{i:03d}.txt"
        path = refs / name
        path.write_text("\n".join(genes) + "\n")
        records.append(record(path))
        references[name] = set(genes)
        lines.append(" ".join(genes))
    unknown = "uncertain_gene"
    (uncertain / "RefOG001.txt").write_text(unknown + "\n")
    records.append(record(uncertain / "RefOG001.txt"))
    lines[0] += " " + unknown
    raw = "\n".join(lines) + "\n"
    tsv = "root_hog\tsource_family\tgenes\n" + "".join(
        f"RootHOG{i:07d}\tFamily{i:07d}\t{','.join(line.split())}\n" for i, line in enumerate(lines))
    old_plain, old_root = tmp_path / "old.txt", tmp_path / "old.tsv"
    old_plain.write_text(raw)
    old_root.write_text(tsv)
    cells = [f"p{p}_c{c}_r{r}" for p, c, r in itertools.product((0, 1), repeat=3)]
    groups = {frozenset(line.split()) for line in lines}
    point = score.score_partition(groups, references, {"RefOG001.txt": {unknown}})
    manifest = dict(references=records,
        predictions={cell: record(old_root if cell.endswith("r1") else old_plain) for cell in cells},
        point_estimates_percent=dict.fromkeys(cells, point))
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    monkeypatch.setattr(score, "ROOT", tmp_path)
    monkeypatch.setattr(score, "FIXED_INPUTS", {"ob_results": ("manifest.json", record(manifest_path)["sha256"])})
    output = tmp_path / "run"
    working = output / "native/orthohmm_working_res"
    phylogeny = output / "native/orthohmm_phylogeny"
    working.mkdir(parents=True)
    phylogeny.mkdir(parents=True)
    plain = working / "orthohmm_edges_clustered.txt"
    root = phylogeny / "orthohmm_root_hogs.tsv"
    plain.write_text(raw)
    root.write_text(tsv)
    return dict(output_root=str(output), genes=141), dict(input_genes=141, orthogroups=70,
        checked_files=[record(plain), record(root)]), manifest, manifest_path


@pytest.mark.parametrize("cell", [f"p{p}_c{c}_r{r}" for p, c, r in itertools.product((0, 1), repeat=3)])
def test_all_factorial_cells_choose_intended_prediction_format(frozen_fixture, cell):
    run, outputs, _, _ = frozen_fixture
    run["cell"] = cell
    evidence = Evidence()
    result = score.score_frozen(run, outputs, evidence)
    assert result["score_fraction"]["f_score"] == 1
    assert result["reference_families"] == 70
    assert result["development_exposed"] is True and result["independent_validation"] is False
    assert result["prediction_format"] == ("root_hogs" if cell.endswith("r1") else "space_separated_groups")
    assert result["canonical_partition_comparison"]["partition_equal"] is True
    assert len(evidence.finish()) == 74


def test_unreviewed_prediction_is_not_used(frozen_fixture):
    run, outputs, _, _ = frozen_fixture
    run["cell"] = "p0_c0_r1"
    outputs["checked_files"].pop()
    with pytest.raises(ValueError, match="not independently output-validated"):
        score.score_frozen(run, outputs, Evidence())


def test_changed_sealed_prediction_refused(frozen_fixture):
    run, outputs, _, _ = frozen_fixture
    run["cell"] = "p0_c0_r0"
    path = Path(outputs["checked_files"][0]["path"])
    path.write_text(path.read_text() + "unknown_gene\n")
    with pytest.raises(ValueError, match="identity changed"):
        score.score_frozen(run, outputs, Evidence())


@pytest.fixture
def joined_scoring(frozen_fixture, monkeypatch):
    run, outputs, _, _ = frozen_fixture
    run.update(index=1, dataset="orthobench", cell="p0_c0_r1", repeat=0)
    root = score.ROOT
    fasta = root / "input.fa"
    genes = [f"g{i:03d}{suffix}" for i in range(1, 71) for suffix in "ab"] + ["uncertain_gene"]
    fasta.write_text("".join(f">{gene}\nAAAA\n" for gene in genes))
    run["inputs"] = [record(fasta)]
    tools = root / "benchmark_tools"
    tools.mkdir()
    original = Path(score.__file__).parent
    for name in ("score_native_orthobench_attempt.py", "score_orthobench_partition.py",
                 "link_factorial_scaling_resources.py", "review_native_factorial_attempt.py"):
        shutil.copyfile(original / name, tools / name)
    baseline = root / "baseline.json"
    baseline.write_text("{}")
    plan_path = root / "plan.json"
    plan_path.write_text(json.dumps(dict(panel_root=str(root / "panel"),
        baseline=record(baseline), helper_sources=[])))
    request_path = root / "request.json"
    request_path.write_text(json.dumps(dict(schema="native_factorial_cost_request_v1", execution_authorized=True,
        job_id=22428, index=1, history=[{}], plan=record(plan_path),
        scheduler_command=str(original / "run_native_factorial_cost.sh"),
        allocation_cwd=str(original.parent))))
    request_ref = record(request_path)
    request = json.loads(request_path.read_text())
    review, _, _, _ = admission_fixture()
    review.update(request=request_ref, plan=record(plan_path), source=record(tools / "review_native_factorial_attempt.py"),
        resources=dict(wall_seconds=1, cpu_seconds=1, peak_memory_bytes=100), timing_disclosure="Shared-host observation.")
    scheduler_path = root / "scheduler.json"
    scheduler_path.write_text("{}")
    review["scheduler"] = record(scheduler_path)
    outputs.update(native_outputs_validated=True, request=request_ref, plan=record(plan_path),
        job_id=22428, index=1, cell=run["cell"], evidence=[])
    outputs_path = root / "outputs.json"
    outputs_path.write_text(json.dumps(outputs))
    review["reviews"] = dict(runtime=record(scheduler_path), resources=record(scheduler_path),
        environment=record(scheduler_path), outputs_or_failure=record(outputs_path))
    review_path = root / "review.json"
    review_path.write_text(json.dumps(review))
    monkeypatch.setattr(score, "validate_plan", lambda plan: [{}, run])
    monkeypatch.setattr(score, "verify_terminal", lambda job: dict(source="fresh_accounting_after_controller_expiry",
        verified=dict(State="COMPLETED", ExitCode="0:0")))
    return request_ref, record(review_path), root / "scoring", outputs_path


def test_joined_score_and_no_overwrite(joined_scoring):
    request_ref, review_ref, destination, _ = joined_scoring
    ref = score.score_attempt(request_ref, review_ref, destination)
    result = json.loads(Path(ref["path"]).read_text())
    assert result["status"] == "terminal_native_orthobench_scored"
    assert result["score_fraction"]["f_score"] == 1
    assert result["native_inference_reexecuted"] is False
    assert result["accuracy_evaluated"] is True
    assert result["publication_ready"] is result["next_identity_authorized"] is False
    assert result["terminal_review"] == review_ref
    with pytest.raises(FileExistsError):
        score.score_attempt(request_ref, review_ref, destination)


def test_changed_actual_scheduler_is_not_scored(joined_scoring, monkeypatch):
    request_ref, review_ref, destination, _ = joined_scoring
    monkeypatch.setattr(score, "verify_terminal", lambda job: dict(source="fresh_accounting_after_controller_expiry",
        verified=dict(State="FAILED", ExitCode="1:0")))
    with pytest.raises(ValueError, match="not successful terminal"):
        score.score_attempt(request_ref, review_ref, destination)
    assert not destination.exists()


def test_live_scheduler_comment_bound(joined_scoring, monkeypatch):
    request_ref, review_ref, destination, _ = joined_scoring
    monkeypatch.setattr(score, "verify_terminal", lambda job: dict(source="live_controller",
        verified=dict(fields=dict(JobState="COMPLETED", ExitCode="0:0", Comment="wrong"))))
    with pytest.raises(ValueError, match="comment differs"):
        score.score_attempt(request_ref, review_ref, destination)
    assert not destination.exists()


def test_changed_reviewed_output_cannot_be_resealed(joined_scoring):
    request_ref, review_ref, destination, output_path = joined_scoring
    output = json.loads(output_path.read_text())
    output["cell"] = "p0_c1_r1"
    output_path.write_text(json.dumps(output))
    review_path = Path(review_ref["path"])
    review = json.loads(review_path.read_text())
    review["reviews"]["outputs_or_failure"] = record(output_path)
    review_path.write_text(json.dumps(review))
    with pytest.raises(ValueError, match="output review identity differs"):
        score.score_attempt(request_ref, record(review_path), destination)
    assert not destination.exists()
