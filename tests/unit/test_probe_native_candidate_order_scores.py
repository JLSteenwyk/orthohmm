"""New native input bindings reuse the original fixed five diagnostic arms."""

import json
from pathlib import Path
import pickle

import numpy as np
import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.probe_native_candidate_order_scores import indexed_seed, load_trusted_hits, prepare, run, validate_species
from benchmark_tools.probe_ob_candidate_order_scores import LABELS


def test_seed_order_preserved(tmp_path):
    path = tmp_path / "seed.txt"
    path.write_text("b a\nc\n")
    assert indexed_seed(path, list("abc")) == [[1, 0], [2]]


@pytest.mark.parametrize("text", ["a b\nb c\n", "a b\n", "a b c d\n", ""])
def test_invalid_seed_membership(tmp_path, text):
    path = tmp_path / "seed.txt"
    path.write_text(text)
    with pytest.raises(ValueError, match="exactly partition"):
        indexed_seed(path, list("abc"))


def test_species_bijection_not_numeric_label_equality():
    validate_species(np.array([0, 1, 0]), np.array([8, 4, 8]))


@pytest.mark.parametrize("fresh", [np.array([8, 8, 8]), np.array([8, 4, 9]),
                                    np.array([8, 4]), np.array([8., 4., 8.])])
def test_species_change_rejected(fresh):
    with pytest.raises(ValueError):
        validate_species(np.array([0, 1, 0]), fresh)


@pytest.mark.parametrize("payload", [{}, [], {"all_gene_ids":["a", "a"], "gene_to_species":{"a":"s"}, "all_hits":{}},
                                     {"all_gene_ids":[], "gene_to_species":{}, "all_hits":{}},
                                     {"all_gene_ids":["a"], "gene_to_species":{"a":"s"}, "all_hits":[]}])
def test_malformed_trusted_cache_rejected(tmp_path, payload):
    path = tmp_path / "fixture.pkl"
    with path.open("wb") as stream:
        pickle.dump(payload, stream)
    with pytest.raises(ValueError):
        load_trusted_hits(record(path))


@pytest.fixture
def prepared(tmp_path, monkeypatch):
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        monkeypatch.setenv(name, "1")

    def write(name, value):
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value) if isinstance(value, (dict, list)) else value)
        return record(path)

    names = list("abcd")
    q, t, s = np.array([0, 2, 0, 3], dtype=np.int32), np.array([2, 0, 3, 0], dtype=np.int32), np.ones(4)
    cache = tmp_path / "cache.pkl"
    with cache.open("wb") as stream:
        pickle.dump({"all_gene_ids": names, "gene_to_species": dict(zip(names, names)),
                     "all_hits": {(names[a], names[b]):float(v) for a, b, v in zip(q, t, s)}}, stream)
    prediction = write("native/prediction.txt", "a c d\nb\n")
    checkpoint = tmp_path / "native/high_sensitivity_checkpoint"
    checkpoint.mkdir()
    arrays = {"gene_to_species.npy": np.array([4, 3, 2, 1], dtype=np.int32),
              "hit_queries.npy": np.r_[q[::-1], np.arange(4)].astype(np.int32),
              "hit_targets.npy": np.r_[t[::-1], np.arange(4)].astype(np.int32),
              "hit_scores.npy": np.r_[s[::-1] + 1e-15, np.ones(4)*2]}
    refs = [write("native/high_sensitivity_checkpoint/gene_names.txt", "a\nb\nc\nd\n")]
    for name, array in arrays.items():
        np.save(checkpoint / name, array, allow_pickle=False)
        refs.append(record(checkpoint / name))
    refs.append(write("native/high_sensitivity_checkpoint/manifest.json",
                      {"schema_version":1, "complete":True, "genes":4, "hits":8}))
    original = write("original.txt", "a c d\nb\n")
    engine = Path(__file__).resolve().parents[2] / "orthohmm/refinement.py"
    prep = {"cache": record(cache), "core_sources": [record(engine)], "candidate_arms": {"p0_c1": {
        "candidate_partition": original, "seed_partition": write("seed.txt", "a\nb\nc d\n"),
        "expansion": {"parameters": {}}}}}
    score = {"schema":"native_factorial_orthobench_score_v1", "status":"terminal_native_orthobench_scored",
             "cell":"p0_c1_r0", "job_id":1, "native_outputs_validated":True,
             "prediction_format":"space_separated_groups", "original_prediction":original,
             "prediction":prediction, "evidence":refs}
    return write, prep, score, tmp_path / "output"


def test_five_arms_end_to_end_and_no_resume(prepared):
    write, prep, score, output = prepared
    plan_ref = prepare(write("preparation.json", prep), write("score.json", score), output)
    result_ref = run(plan_ref)
    result = json.loads(Path(result_ref["path"]).read_text())
    assert result["status"] == "controls_reproduced"
    assert result["controls"] == {"historical":True, "fresh_full":True}
    assert [r["label"] for r in result["rows"]] == list(LABELS)
    assert result["removed_self_hits"] == 4
    assert len(result["comparisons"]) == 10
    assert not result["accuracy_scored"]
    with pytest.raises(FileExistsError, match="do not retry"):
        run(plan_ref)


def test_wrong_native_control_retains_all_arms_and_failure(prepared):
    write, prep, score, output = prepared
    score["prediction"] = write("native/prediction.txt", "a b c d\n")
    plan_ref = prepare(write("preparation.json", prep), write("score.json", score), output)
    with pytest.raises(ValueError, match="control mismatch"):
        run(plan_ref)
    failure = json.loads((output / "failure.json").read_text())
    assert failure["completed_labels"] == list(LABELS)
    assert failure["retry"] is False
    assert json.loads((output / "report.json").read_text())["status"] == "control_mismatch"


def test_changed_gene_indexing_fails_before_any_arm(prepared):
    write, prep, score, output = prepared
    changed = write("native/high_sensitivity_checkpoint/gene_names.txt", "b\na\nc\nd\n")
    score["evidence"] = [changed if r["path"] == changed["path"] else r for r in score["evidence"]]
    plan_ref = prepare(write("preparation.json", prep), write("score.json", score), output)
    with pytest.raises(ValueError, match="indexing differs"):
        run(plan_ref)
    assert json.loads((output / "failure.json").read_text())["completed_labels"] == []


@pytest.mark.parametrize("change", ["missing_checkpoint", "unsupported_cell", "wrong_original"])
def test_prepare_rejects_incompatible_boundaries(prepared, change):
    write, prep, score, output = prepared
    if change == "missing_checkpoint":
        score["evidence"].pop()
    elif change == "unsupported_cell":
        score["cell"] = "p0_c1_r1"
    else:
        score["original_prediction"] = record(__file__)
    with pytest.raises(ValueError):
        prepare(write("preparation.json", prep), write("score.json", score), output)
    assert not output.exists()


def test_runtime_change_rejected_before_started(prepared, monkeypatch):
    write, prep, score, output = prepared
    plan_ref = prepare(write("preparation.json", prep), write("score.json", score), output)
    monkeypatch.setenv("OMP_NUM_THREADS", "2")
    with pytest.raises(ValueError, match="thread settings changed"):
        run(plan_ref)
    assert not (output / "started.json").exists()


def test_changed_bound_input_rejected_before_started(prepared):
    write, prep, score, output = prepared
    plan_ref = prepare(write("preparation.json", prep), write("score.json", score), output)
    Path(prep["candidate_arms"]["p0_c1"]["seed_partition"]["path"]).write_text("a b c d\n")
    with pytest.raises(ValueError, match="identity changed"):
        run(plan_ref)
    assert not (output / "started.json").exists()
