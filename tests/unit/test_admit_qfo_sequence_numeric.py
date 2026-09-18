from pathlib import Path

import pytest

from benchmark_tools import admit_qfo_sequence_numeric as module


def fixture():
    root = Path("/numeric")
    variants = {}
    for label, cap in (("all_hits", None), ("top100", 100)):
        path = root / label / "orthohmm_working_res/high_sensitivity_checkpoint"
        manifest = {"path": str(path / "manifest.json"), "bytes": 1, "sha256": "fixture"}
        variants[label] = {"checkpoint": str(path), "cap": cap, "manifest": manifest,
                           "audit": {"status": "numeric_checkpoint_verified", "manifest": manifest}}
    source = {"path": "/source", "bytes": 1, "sha256": "fixture"}
    report = {"status": "corrected_qfo_sequence_numeric_checkpoints_verified", "source": source,
        "job_id": "123", "numeric_validated": True, "accuracy_evaluated": False, "publication_ready": False,
        "genes": 984137, "proteomes": 78, "hits": 100, "variants": variants}
    return report, {"JobIDRaw": "123"}, source, root


@pytest.mark.parametrize("problem", [None, "status", "job", "source", "genes", "species", "hits",
    "missing", "path", "cap", "float_cap", "audit", "accuracy"])
def test_conversion_identity(problem):
    report, scheduler, source, root = fixture()
    if problem == "status":
        report["status"] = "failed"
    elif problem == "job":
        report["job_id"] = "999"
    elif problem == "source":
        report["source"] = {}
    elif problem == "genes":
        report["genes"] = 976504
    elif problem == "species":
        report["proteomes"] = 12
    elif problem == "hits":
        report["hits"] = True
    elif problem == "missing":
        report["variants"].pop("top100")
    elif problem == "path":
        report["variants"]["all_hits"]["checkpoint"] = "/historical"
    elif problem == "cap":
        report["variants"]["all_hits"]["cap"] = 100
    elif problem == "float_cap":
        report["variants"]["top100"]["cap"] = 100.0
    elif problem == "audit":
        report["variants"]["top100"]["audit"]["manifest"] = {}
    elif problem == "accuracy":
        report["accuracy_evaluated"] = True
    if problem:
        with pytest.raises(ValueError):
            module.validate_conversion(report, scheduler, source, root)
    else:
        module.validate_conversion(report, scheduler, source, root)


@pytest.mark.parametrize("state,cpus,memory", [("COMPLETED", "2", "192G"), ("RUNNING", "2", "192G"),
    ("FAILED", "2", "192G"), ("COMPLETED", "32", "192G"), ("COMPLETED", "2", "64G")])
def test_terminal_allocation(monkeypatch, state, cpus, memory):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS|ReqMem\n" + f"123|{state}|0:0|01:00|bizon|{cpus}|{memory}\n")
    if (state, cpus, memory) == ("COMPLETED", "2", "192G"):
        assert module.completed("123", 2, "192G")["JobIDRaw"] == "123"
    else:
        with pytest.raises(ValueError):
            module.completed("123", 2, "192G")


def test_pending_conversion_does_not_read_files(tmp_path, monkeypatch):
    def pending(*args):
        raise ValueError("Pending")
    monkeypatch.setattr(module, "completed", pending)
    monkeypatch.setattr(module, "record", lambda *a: pytest.fail("Read pending conversion"))
    with pytest.raises(ValueError, match="Pending"):
        module.admit(tmp_path, "123", tmp_path / "audit")
    assert not (tmp_path / "audit").exists()
