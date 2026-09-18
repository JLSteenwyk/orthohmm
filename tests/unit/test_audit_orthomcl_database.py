import json
import os
from types import SimpleNamespace

import pytest

from benchmark_tools import audit_orthomcl_database as module


def fixture(tmp_path):
    fasta, dump = tmp_path / "all.fa", tmp_path / "database.fa"
    fasta.write_text(">sp|A|A_ONE\nACDE\n>B description\nMUA\n")
    dump.write_text(">gnl|BL_ORD_ID|0 sp|A|A_ONE\nACDE\n>gnl|BL_ORD_ID|1 B description\nMUA\n")
    return fasta, dump


def test_exact_parity(tmp_path):
    result = module.compare_dump(*fixture(tmp_path))
    assert result["exact_sequence_parity"] is True
    assert result["input_sequences"] == result["exact_sequence_matches"] == 2
    assert result["input_residues"] == result["database_residues"] == 7


@pytest.mark.parametrize("sequence", ["MXA", "mua", "MU", "MUAA"])
def test_residue_and_length_changes_reported_without_normalization(tmp_path, sequence):
    fasta, dump = fixture(tmp_path)
    dump.write_text(dump.read_text().replace("MUA", sequence))
    result = module.compare_dump(fasta, dump)
    assert result["exact_sequence_parity"] is False
    assert result["sequence_difference_count"] == 1
    assert result["differences"][0]["id"] == "B"


@pytest.mark.parametrize("replacement", ["gnl|BL_ORD_ID|5", "lcl|0", "gnl|BL_ORD_ID|1"])
def test_changed_ordinal_header_rejected(tmp_path, replacement):
    fasta, dump = fixture(tmp_path)
    dump.write_text(dump.read_text().replace("gnl|BL_ORD_ID|0", replacement))
    with pytest.raises(ValueError, match="header/order"):
        module.compare_dump(fasta, dump)


@pytest.mark.parametrize("kind", ["missing", "extra", "duplicate", "empty"])
def test_cardinality_and_duplicate_rejected(tmp_path, kind):
    fasta, dump = fixture(tmp_path)
    if kind == "missing":
        dump.write_text(">gnl|BL_ORD_ID|0 sp|A|A_ONE\nACDE\n")
    elif kind == "extra":
        dump.write_text(dump.read_text() + ">gnl|BL_ORD_ID|2 C\nMA\n")
    elif kind == "duplicate":
        fasta.write_text(">A\nMA\n>A\nMA\n")
        dump.write_text(">gnl|BL_ORD_ID|0 A\nMA\n>gnl|BL_ORD_ID|1 A\nMA\n")
    else:
        fasta.write_text("")
        dump.write_text("")
    with pytest.raises(ValueError):
        module.compare_dump(fasta, dump)


def make_database(fasta):
    for suffix in ("phr", "pin", "psq"):
        fasta.with_name(fasta.name + "." + suffix).write_text("fixture")


@pytest.mark.parametrize("problem", ["missing", "empty", "alias", "symlink"])
def test_bad_database_inventory(tmp_path, problem):
    fasta, _ = fixture(tmp_path)
    make_database(fasta)
    path = tmp_path / "all.fa.pin"
    if problem == "missing":
        path.unlink()
    elif problem == "empty":
        path.write_text("")
    elif problem == "alias":
        (tmp_path / "all.fa.pal").write_text("alias")
    else:
        path.unlink()
        path.symlink_to(fasta)
    with pytest.raises(ValueError, match="three nonempty"):
        module.database_files(fasta)


@pytest.mark.parametrize("problem", [None, "exit", "diagnostics", "mutation", "difference"])
def test_audit_preserves_report_and_checks_inputs(tmp_path, monkeypatch, problem):
    fasta, dump = fixture(tmp_path)
    make_database(fasta)
    tool = tmp_path / "fastacmd"
    tool.write_text("fixture")
    runtime_path = tmp_path / "runtime.json"
    runtime_path.write_text("fixture")
    runtime = {"records": [{**module.record(tool), "kind": "file"}]}
    monkeypatch.setattr(module, "FASTACMD", tool)
    monkeypatch.setattr(module, "read_frozen", lambda *a: runtime)
    monkeypatch.setattr(module, "verify_tree", lambda *a: None)
    monkeypatch.setattr(module.Path, "home", lambda: tmp_path)
    def run(argv, **kwargs):
        text = dump.read_text().replace("MUA", "MXA") if problem == "difference" else dump.read_text()
        kwargs["stdout"].write(text.encode())
        if problem == "diagnostics":
            kwargs["stderr"].write(b"unexpected warning")
        if problem == "mutation":
            fasta.write_text(">changed\nMA\n")
        return SimpleNamespace(returncode=1 if problem == "exit" else 0)
    monkeypatch.setattr(module.subprocess, "run", run)
    output = tmp_path / "audit"
    if problem in ("exit", "diagnostics", "mutation"):
        with pytest.raises(ValueError):
            module.audit(fasta, runtime_path, output)
    else:
        result = module.audit(fasta, runtime_path, output)
        assert result["content"]["exact_sequence_parity"] == (problem is None)
    report = json.loads((output / "report.json").read_text())
    assert report["database_admitted"] is False and report["search_admitted"] is False
    assert len(report["outputs"]) == 2
    assert (report["status"] == "failed") == (problem in ("exit", "diagnostics", "mutation"))
    with pytest.raises(FileExistsError):
        module.audit(fasta, runtime_path, output)


@pytest.mark.skipif(os.environ.get("ORTHOHMM_LEGACY_BLAST_SMOKE") != "1",
                    reason="Opt-in installed legacy database round trip")
def test_installed_formatdb_normalization_is_detected(tmp_path):
    fasta = tmp_path / "all.fa"
    fasta.write_text(">sp|PROBE1|CANONICAL\nACDEFGHIKLMNPQRSTVWY\n"
                     ">sp|PROBE2|AMBIGUOUS\nMABZXU*\n>sp|PROBE3|LOWERCASE\nmacdefghik\n")
    module.subprocess.run([str(module.FASTACMD.with_name("formatdb")), "-i", str(fasta), "-p", "t"],
                          cwd=tmp_path, env=module.environment(), check=True, capture_output=True, timeout=30)
    dump = tmp_path / "extracted.fasta"
    with dump.open("wb") as stream:
        module.subprocess.run([str(module.FASTACMD), "-d", str(fasta), "-p", "T", "-D", "1"],
                              cwd=tmp_path, env=module.environment(), stdout=stream,
                              stderr=module.subprocess.PIPE, check=True, timeout=30)
    content = module.compare_dump(fasta, dump)
    assert content["sequence_difference_count"] == 2
    assert content["input_residues"] == 37 and content["database_residues"] == 36
    assert content["exact_sequence_parity"] is False
    assert "MABZXU\n" in dump.read_text()
    assert "MACDEFGHIK\n" in dump.read_text()
