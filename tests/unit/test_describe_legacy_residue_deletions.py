import copy
import io
import json
import os
import subprocess

import pytest
from Bio.Blast import NCBIXML

from benchmark_tools import describe_legacy_residue_deletions as module
from benchmark_tools.audit_orthomcl_database import FASTACMD, environment


def fixture(tmp_path, source="MAOCO", native="MAC"):
    fasta, dump = tmp_path / "all.fa", tmp_path / "database.fasta"
    fasta.write_text(f">A\n{source}\n>B\nMXA\n")
    dump.write_text(f">gnl|BL_ORD_ID|0 A\n{native}\n>gnl|BL_ORD_ID|1 B\nMXA\n")
    return fasta, dump, module.compare_dump(fasta, dump)


def test_positions_and_no_admission(tmp_path):
    result = module.describe(*fixture(tmp_path))
    assert result["transformations"][0]["positions_one_based"] == [3, 5]
    assert result["content"]["exact_sequence_matches"] == 1
    for key in ("database_admitted", "search_admitted", "publication_ready", "exact_sequence_parity"):
        assert result[key] is False


@pytest.mark.parametrize("source,native", [("MAC", "MA"), ("MAOC", "MAXC"),
                                          ("MAOC", "MC"), ("MAOC", "MAOC")])
def test_other_changes_and_exact_parity_rejected(tmp_path, source, native):
    with pytest.raises(ValueError):
        module.describe(*fixture(tmp_path, source, native))


@pytest.mark.parametrize("key", ["input_residues", "exact_sequence_matches", "differences"])
def test_changed_reviewed_content_rejected(tmp_path, key):
    fasta, dump, content = fixture(tmp_path)
    changed = copy.deepcopy(content)
    changed[key] = [] if key == "differences" else -1
    with pytest.raises(ValueError, match="reviewed evidence"):
        module.describe(fasta, dump, changed)


def test_unreviewed_extra_change_rejected(tmp_path):
    fasta, dump, content = fixture(tmp_path)
    dump.write_text(dump.read_text().replace("MXA", "MA"))
    with pytest.raises(ValueError, match="reviewed evidence"):
        module.describe(fasta, dump, content)


@pytest.mark.parametrize("problem", [None, "digest", "source", "dump", "status", "missing_identity"])
def test_audit_checks_bound_evidence(tmp_path, problem):
    fasta, dump, content = fixture(tmp_path)
    report = dict(status="database_sequence_differences_require_review",
                  command=["fastacmd", "-d", str(fasta)], content=content,
                  checked_records=[module.record(fasta)], outputs=[module.record(dump)])
    if problem == "status":
        report["status"] = "running"
    if problem == "missing_identity":
        report["outputs"] = []
    path = tmp_path / "report.json"
    path.write_text(json.dumps(report))
    sha = module.record(path)["sha256"]
    if problem == "digest":
        sha = "0" * 64
    if problem == "source":
        fasta.write_text(">changed\nAA\n")
    if problem == "dump":
        dump.write_text(">changed\nAA\n")
    output = tmp_path / "description.json"
    if problem:
        with pytest.raises(ValueError):
            module.audit(path, sha, output)
        assert not output.exists()
    else:
        result = module.audit(path, sha, output)
        assert json.loads(output.read_text()) == result
        with pytest.raises(FileExistsError):
            module.audit(path, sha, output)


@pytest.mark.skipif(os.environ.get("ORTHOHMM_LEGACY_BLAST_SMOKE") != "1",
                   reason="Opt-in installed native query and database parser probe")
def test_native_query_and_database_delete_o(tmp_path):
    with_o = "ACDEFGHIKLMNOPQRSTVWY"
    fasta = tmp_path / "input.fa"
    fasta.write_text(f">with_o\n{with_o}\n>without_o\n{with_o.replace('O', '')}\n"
                     f">with_x\n{with_o.replace('O', 'X')}\n")
    done = subprocess.run([str(FASTACMD.with_name("formatdb")), "-i", str(fasta), "-p", "T"],
                          cwd=tmp_path, env=environment(), capture_output=True, check=True, timeout=30)
    assert b"1 illegal character was removed" in done.stderr and b"1 O" in done.stderr
    dumped = subprocess.run([str(FASTACMD), "-d", str(fasta), "-p", "T", "-D", "1"],
                            cwd=tmp_path, env=environment(), capture_output=True, check=True, timeout=30)
    dump = tmp_path / "database.fasta"
    dump.write_bytes(dumped.stdout)
    content = module.compare_dump(fasta, dump)
    assert content["sequence_difference_count"] == 1
    assert module.describe(fasta, dump, content)["transformations"][0]["positions_one_based"] == [13]
    search = subprocess.run([str(FASTACMD.with_name("blastall")), "-p", "blastp", "-i", str(fasta),
                             "-d", str(fasta), "-F", "F", "-e", "10", "-m", "7"],
                            cwd=tmp_path, env=environment(), capture_output=True, check=True, timeout=30)
    queries = list(NCBIXML.parse(io.StringIO(search.stdout.decode())))
    assert [q.query for q in queries] == ["with_o", "without_o", "with_x"]
    assert [q.query_length for q in queries] == [len(with_o) - 1, len(with_o) - 1, len(with_o)]
    assert any(h.query == with_o.replace("O", "") for a in queries[0].alignments for h in a.hsps)
