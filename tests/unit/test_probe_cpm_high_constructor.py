import copy
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from benchmark_tools.probe_cpm_high_constructor import require_audit, require_result, run


def test_actual_audit_and_changed_target():
    root = Path(__file__).resolve().parents[2]
    report = json.loads((root / "benchmark_tools/results/qfo_cpm_high_failed_payload_audit_20260923.json").read_text())
    require_audit(report)
    report["observed"]["edges"] -= 1
    with pytest.raises(ValueError):
        require_audit(report)


@pytest.mark.parametrize("change", [None, "weight", "before", "after", "optimizer", "format"])
def test_result_gate(change):
    differences = {k: {"different_edges": 0} for k in
                   ("native_vs_saved", "constructor_vs_saved", "native_vs_constructor")}
    report = dict(status="direct_construction_observed", mode="minimal_imports", edge_format="python_pairs",
                  optimizer_called=False, accuracy_evaluated=False, saved={"digest": "a"},
                  before_weights=dict(differences=copy.deepcopy(differences)),
                  after_weights=dict(fingerprint={"digest": "a"}, differences=copy.deepcopy(differences)))
    if change == "weight":
        report["after_weights"]["fingerprint"]["digest"] = "b"
    elif change in ("before", "after"):
        report[change + "_weights"]["differences"]["native_vs_saved"]["different_edges"] = 1
    elif change == "optimizer":
        report["optimizer_called"] = True
    elif change == "format":
        report["edge_format"] = "numpy"
    if change:
        with pytest.raises(ValueError):
            require_result(report)
    else:
        require_result(report)


def test_no_overwrite_or_dangling_symlink(tmp_path):
    with pytest.raises(FileExistsError):
        run(tmp_path, tmp_path)
    link = tmp_path / "link"
    link.symlink_to(tmp_path / "missing")
    with pytest.raises(FileExistsError):
        run(tmp_path, link)


@pytest.mark.parametrize("mode", ["minimal_imports", "frozen_imports"])
def test_fresh_native_worker_on_small_graph(tmp_path, mode):
    root = Path(__file__).resolve().parents[2]
    payload = tmp_path / "payload"
    payload.mkdir()
    (payload / "gene_names.txt").write_text("a\nb\nc\nd\n")
    for name, values in (("sources", np.array([2, 0, 1], dtype=np.int32)),
                         ("targets", np.array([1, 2, 0], dtype=np.int32)),
                         ("weights", np.array([.3, 1.2, 2.5], dtype=np.float64))):
        np.save(payload / (name + ".npy"), values)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONNOUSERSITE="1")
    subprocess.run([sys.executable, "-B", str(root / "benchmark_tools/probe_cpm_high_constructor.py"),
                    "--root", str(root), "--output", str(tmp_path), "--worker-payload", str(payload), "--mode", mode],
                   check=True, env=env, capture_output=True, text=True, timeout=60)
    result = json.loads((payload / "result.json").read_text())
    require_result(result, mode)
    assert result["saved"]["vertices"] == 4
    assert result["saved"]["edges"] == 3
    before = json.loads((payload / "worker_before.json").read_text())
    imported = {k.split(".")[0] for k in before["modules"]} & {"orthohmm", "leidenalg"}
    assert imported == (set() if mode == "minimal_imports" else {"orthohmm", "leidenalg"})


def test_unknown_mode_rejected(tmp_path):
    with pytest.raises(ValueError, match="Unknown import mode"):
        run(tmp_path, tmp_path / "out", "other")
