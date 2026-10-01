import subprocess
import sys

import pytest

from benchmark_tools import readback_qfo_private_phylogeny_control as module


def accounting():
    return ("JobID|JobIDRaw|State|ExitCode|Elapsed|AllocCPUS|ReqMem|NodeList\n"
        "22387|22387|COMPLETED|0:0|00:10:28|32|192G|bizon\n"
        "22387.batch|22387.batch|COMPLETED|0:0|00:10:28|32||bizon\n"
        "22389|22389|COMPLETED|0:0|00:04:04|2|64G|bizon\n"
        "22389.batch|22389.batch|COMPLETED|0:0|00:04:04|2||bizon\n")


@pytest.mark.parametrize("old,new", [("COMPLETED", "RUNNING"), ("COMPLETED", "FAILED"),
    ("0:0", "1:0"), ("32|192G", "2|64G"), ("2|64G", "32|192G"), ("bizon", "dgx")])
def test_completion_gate(old, new):
    with pytest.raises(ValueError):
        module.completed(accounting().replace(old, new))


@pytest.mark.parametrize("problem", ["missing_parent", "missing_step", "duplicate_parent", "failed_step"])
def test_partial_completion(problem):
    text = accounting()
    if problem == "missing_parent":
        text = "\n".join(line for line in text.splitlines() if not line.startswith("22389|"))
    elif problem == "missing_step":
        text = "\n".join(line for line in text.splitlines() if not line.startswith("22389.batch|"))
    elif problem == "duplicate_parent":
        text += text.splitlines()[1] + "\n"
    else:
        text = text.replace("22389.batch|22389.batch|COMPLETED", "22389.batch|22389.batch|FAILED")
    with pytest.raises(ValueError):
        module.completed(text)


def test_valid_completion():
    assert [row["JobID"] for row in module.completed(accounting())] == ["22387", "22389"]


def test_hash_record(tmp_path):
    path = tmp_path / "file"
    path.write_bytes(b"abc")
    assert module.record(path) == dict(path=str(path), bytes=3,
        sha256="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")


def test_readback_has_no_scientific_imports():
    code = ("import sys; import benchmark_tools.readback_qfo_private_phylogeny_control; "
            "assert not any(n == 'orthohmm' or n.startswith('orthohmm.') or n.startswith('numpy') for n in sys.modules)")
    subprocess.run([sys.executable, "-S", "-B", "-c", code], check=True)
