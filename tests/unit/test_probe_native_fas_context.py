from copy import deepcopy
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from benchmark_tools import probe_native_fas_context as current


def rows():
    data = []
    for context in current.contexts():
        raw = {"_".join(pair): (["NA", "NA"] if pair == current.FOCAL[3] else
                               ["1", "1"] if pair == current.FOCAL[0] else ["0.6", "0.8"])
               for pair in context["pairs"]}
        data.append(dict(context={**context, "pairs": [list(p) for p in context["pairs"]]},
                         raw=raw, normalized=current.canonical_raw(raw, context["pairs"]),
                         terminal=dict(returncode=0), error=None))
    return data


def test_full_plan_covers_numeric_na_companions_order_workers_and_hashes():
    plan = current.contexts()
    assert len(plan) == 12 and sum(len(c["pairs"]) for c in plan) == 40
    assert {c["workers"] for c in plan} == {1, 2} and {c["hash_seed"] for c in plan} == {"0", "1"}
    result = current.compare(rows())
    assert result["complete"] and result["all_tested_contexts_invariant"] and result["pair_evaluations"] == 40


def test_fixture_annotations_do_not_share_mutable_gene_state():
    annotations, owners = current.fixtures()
    assert set(owners) == {gene for doc in annotations.values() for gene in doc["feature"]}
    a = annotations["T1"]["feature"]["A"]
    d = deepcopy(annotations["T2"]["feature"]["D"])
    a["Pfam"].clear()
    assert annotations["T2"]["feature"]["D"] == d
    assert annotations["T1"]["feature"]["C"] != annotations["T2"]["feature"]["F"]


@pytest.mark.parametrize("raw", [{}, {"A_D": ["0.2"]}, {"A_D": ["NaN", ".2"]}, {"A_D": ["1.2", ".2"]}])
def test_incomplete_or_invalid_native_json_refused(raw):
    with pytest.raises(ValueError):
        current.canonical_raw(raw, [current.FOCAL[0]])


def test_mixed_na_is_omitted_not_half_direction_score():
    result = current.canonical_raw({"A_D": ["NA", "0.8"]}, [current.FOCAL[0]])
    assert result["A_D"]["value"] is None and not result["A_D"]["returned"]


def test_context_sensitive_scores_are_reported_not_accepted():
    data = rows()
    data[4]["raw"]["B_E"] = ["0.5", "0.5"]
    data[4]["normalized"] = current.canonical_raw(data[4]["raw"], current.contexts()[4]["pairs"])
    result = current.compare(data)
    assert result["complete"] and not result["all_tested_contexts_invariant"]
    assert result["variations"][0]["pair"] == ["B", "E"]


def test_failed_or_incomplete_contexts_are_not_admitted():
    data = rows()
    data[-1]["terminal"]["returncode"] = 1
    assert not current.compare(data)["all_tested_contexts_invariant"]
    assert not current.compare(data[:-1])["all_tested_contexts_invariant"]


@pytest.mark.parametrize("change", ["plan", "raw_binding", "identical", "nonidentical", "priority", "na"])
def test_changed_plan_bindings_or_controls_refused(change):
    data = rows()
    if change == "plan": data[0]["context"]["workers"] = 2
    elif change == "raw_binding": data[0]["normalized"]["A_D"]["value"] = .3
    else:
        index, key, values = {"identical": (0, "A_D", [".3", ".3"]),
                              "nonidentical": (1, "B_E", ["1", "1"]),
                              "priority": (2, "C_F", ["NA", "NA"]),
                              "na": (3, "G_H", [".2", ".2"])}[change]
        data[index]["raw"][key] = values
        data[index]["normalized"] = current.canonical_raw(data[index]["raw"], current.contexts()[index]["pairs"])
    with pytest.raises(ValueError):
        current.compare(data)


def test_native_io_capture_uses_real_helper_return_without_replacement(tmp_path, monkeypatch):
    context = current.contexts()[0]
    raw = {"A_D": ["1", "1"]}
    out = tmp_path / "fas"
    out.mkdir()
    def child(argv, seed):
        (out / "computed_results.json").write_text(json.dumps(raw))
        return current.subprocess.CompletedProcess(argv, 0, "native output", ""), False
    monkeypatch.setattr(current, "run_child", child)
    def helper(pairs, owners, annotations, workers):
        scorer.subprocess.run(["fas.runMultiTaxa", "-o", out])
        return {("A", "D"): 1.}
    scorer = SimpleNamespace(subprocess=current.subprocess, compute_fas_scores_for_pairs=helper)
    result = current.execute_context(scorer, context, tmp_path, {})
    assert result["raw"] == raw and result["normalized"]["A_D"]["value"] == 1
    assert result["terminal"]["stdout"] == "native output" and result["error"] is None


def test_existing_output_refused_before_native_imports_or_affinity(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["probe", "--protocol", "missing", "--output", str(tmp_path)])
    with pytest.raises(FileExistsError):
        current.main()


@pytest.mark.parametrize("payload", ["not JSON", '{"A_D": [NaN, 0.5]}', '{}', '{"A_D": ["NaN", ".5"]}', None])
def test_invalid_native_output_preserves_terminal_and_original_text(tmp_path, monkeypatch, payload):
    def child(argv, seed):
        if payload is not None:
            (tmp_path / "computed_results.json").write_text(payload)
        return current.subprocess.CompletedProcess(argv, 0, "retained stdout", "retained stderr"), False
    monkeypatch.setattr(current, "run_child", child)
    def helper(pairs, owners, annotations, workers):
        scorer.subprocess.run(["fas.runMultiTaxa", "-o", tmp_path])
        return {}
    scorer = SimpleNamespace(subprocess=current.subprocess, compute_fas_scores_for_pairs=helper)
    row = current.execute_context(scorer, current.contexts()[0], tmp_path, {})
    assert row["error"] is not None and row["raw_text"] == payload
    assert row["terminal"]["returncode"] == 0 and row["terminal"]["stdout"] == "retained stdout"
    assert row["terminal"]["stderr"] == "retained stderr"
    assert not current.compare([row])["all_tested_contexts_invariant"]
    json.dumps(row, allow_nan=False)


def test_helper_failure_before_child_is_retained(tmp_path):
    def helper(*args):
        raise RuntimeError("before child launch")
    scorer = SimpleNamespace(subprocess=current.subprocess, compute_fas_scores_for_pairs=helper)
    row = current.execute_context(scorer, current.contexts()[0], tmp_path, {})
    assert row["terminal"]["returncode"] is None and row["terminal"]["command"] is None
    assert row["error"]["message"] == "before child launch"
    assert not current.compare([row])["all_tested_contexts_invariant"]


def test_child_launch_failure_preserves_attempted_command(tmp_path, monkeypatch):
    def child(*args):
        raise FileNotFoundError("child executable missing")
    monkeypatch.setattr(current, "run_child", child)
    def helper(*args):
        scorer.subprocess.run(["fas.runMultiTaxa", "-o", tmp_path])
    scorer = SimpleNamespace(subprocess=current.subprocess, compute_fas_scores_for_pairs=helper)
    row = current.execute_context(scorer, current.contexts()[0], tmp_path, {})
    assert row["terminal"]["command"][0] == "fas.runMultiTaxa"
    assert row["error"]["type"] == "FileNotFoundError"


def test_native_loader_mismatch_is_retained_not_thrown(tmp_path, monkeypatch):
    def child(argv, seed):
        (tmp_path / "computed_results.json").write_text('{"A_D": ["1", "1"]}')
        return current.subprocess.CompletedProcess(argv, 0, "", ""), False
    monkeypatch.setattr(current, "run_child", child)
    def helper(*args):
        scorer.subprocess.run(["fas.runMultiTaxa", "-o", tmp_path])
        return {("A", "D"): .5}
    scorer = SimpleNamespace(subprocess=current.subprocess, compute_fas_scores_for_pairs=helper)
    row = current.execute_context(scorer, current.contexts()[0], tmp_path, {})
    assert row["raw"] == {"A_D": ["1", "1"]}
    assert row["error"]["stage"] == "native_output_validation"
    assert "Native loader differs" in row["error"]["message"]


def test_timeout_kills_and_reaps_only_own_new_process_group(monkeypatch):
    calls = []
    class Child:
        pid = 12345
        returncode = -9
        def communicate(self, timeout=None):
            calls.append(timeout)
            if timeout is not None:
                raise current.subprocess.TimeoutExpired(["native"], timeout)
            return "partial stdout", "partial stderr"
    def popen(command, **kwargs):
        assert kwargs["start_new_session"] and kwargs["env"]["PYTHONHASHSEED"] == "1"
        return Child()
    monkeypatch.setattr(current.subprocess, "Popen", popen)
    killed = []
    monkeypatch.setattr(current.os, "killpg", lambda pid, signal: killed.append((pid, signal)))
    terminal, timed_out = current.run_child(["native"], "1")
    assert timed_out and terminal.returncode == -9
    assert terminal.stdout == "partial stdout" and terminal.stderr == "partial stderr"
    assert calls == [120, None] and killed == [(12345, current.signal.SIGKILL)]


def test_protocol_is_prospectively_bound():
    root = Path(__file__).resolve().parents[2]
    assert current.fingerprint(root / "benchmark_tools/results/NATIVE_FAS_SAMPLING_DESIGN_PROTOCOL_20261010.md")["sha256"] == current.PROTOCOL_SHA
