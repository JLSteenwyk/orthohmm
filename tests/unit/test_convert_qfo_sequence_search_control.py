from pathlib import Path
import sqlite3

import pytest

from benchmark_tools import convert_qfo_sequence_search_control as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from orthohmm.accuracy import load_accuracy_checkpoint


@pytest.mark.parametrize("problem", [None, "length", "bool", "species", "extra", "missing", "duplicate", "scope"])
def test_metadata_bound_to_fasta(tmp_path, problem):
    fasta = tmp_path / "s.fa"
    fasta.write_text(">a\nAAAA\n>b\nAAA\n" + (">a\nAAAA\n" if problem == "duplicate" else ""))
    metadata = {"a": {"length": 4, "species": "s.fa"}, "b": {"length": 3, "species": "s.fa"}}
    if problem == "length":
        metadata["a"]["length"] = 5
    elif problem == "bool":
        metadata["a"]["length"] = True
    elif problem == "species":
        metadata["a"]["species"] = "wrong"
    elif problem == "extra":
        metadata["c"] = {"length": 3, "species": "s.fa"}
    elif problem == "missing":
        metadata.pop("b")
    args = ([record(fasta)], metadata, 3 if problem == "scope" else 2, 1)
    if problem:
        with pytest.raises(ValueError):
            module.validate_metadata(*args)
    else:
        module.validate_metadata(*args)


def test_numeric_variants_keep_all_hits_and_cap_only_diagnostic(tmp_path):
    metadata = {"q": {"length": 20, "species": "a.fa"},
                **{f"t{i:03d}": {"length": 20, "species": "b.fa"} for i in range(101)}}
    own, other = tmp_path / "own.tsv", tmp_path / "other.tsv"
    own.write_text("q\tq\t20\t20\t40\t20\t1e-10\n")
    other.write_text("".join(f"q\tt{i:03d}\t20\t20\t40\t20\t1e-10\n" for i in reversed(range(101))))
    plan = {"searches": [{"target_fasta": {"path": "/a.fa"}, "output": str(own)},
                         {"target_fasta": {"path": "/b.fa"}, "output": str(other)}]}
    output = tmp_path / "numeric"
    output.mkdir()
    variants, count = module.convert_hits(plan, metadata, output)
    assert count == 102
    assert variants["all_hits"]["audit"]["summary"]["hits"] == 102
    assert variants["top100"]["audit"]["summary"]["hits"] == 101
    assert all(v["audit"]["summary"]["self_hits"] == 1 for v in variants.values())
    names, species, queries, targets, scores = load_accuracy_checkpoint(variants["top100"]["checkpoint"])
    assert set(scores) == {2.}
    assert "t100" not in {names[i] for i in targets}
    assert set(queries) == {names.index("q")}
    assert set(species) == {0, 1}


@pytest.mark.parametrize("line,error", [("q\tq\t20\t20\t40\t20\t1e-10\n" * 2, sqlite3.IntegrityError),
    ("q\tq\t21\t20\t40\t20\t1e-10\n", ValueError), ("", ValueError)])
def test_invalid_rows_cannot_produce_variants(tmp_path, line, error):
    hits = tmp_path / "hits.tsv"
    hits.write_text(line)
    output = tmp_path / "numeric"
    output.mkdir()
    with pytest.raises(error):
        module.convert_hits({"searches": [{"target_fasta": {"path": "/a.fa"}, "output": str(hits)}]},
                            {"q": {"length": 20, "species": "a.fa"}}, output)
    assert not (output / "all_hits").exists()


def test_pending_admission_prevents_partial_reads(tmp_path, monkeypatch):
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS|ReqMem\n456|PENDING|0:0|00:00|bizon|2|64G\n")
    monkeypatch.setattr(module, "read_frozen", lambda *a: pytest.fail("Read pending evidence"))
    with pytest.raises(ValueError):
        module.convert(tmp_path, tmp_path / "absent", "0" * 64, "456", tmp_path / "numeric")
    assert not (tmp_path / "numeric").exists()
