import json

import pytest

from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.prepare_ygob_validation import prepare


def data(tmp_path):
    candidate = tmp_path / "candidate"
    candidate.mkdir()
    fields = ["---"] * 33
    fields[0], fields[11], fields[13] = "a", "exposed", "b"
    (candidate / "Pillars.tab").write_text("\t".join(fields) + "\n")
    (candidate / "AA.fsa").write_text(">a {ON}\nAAAA*\n>b {ON}\nC*CC\n>exposed {ON}\nAAAA\n>off {OFF}\nCCCC\n")
    (candidate / "README").write_text("test fixture")
    audit = tmp_path / "audit.json"
    audit.write_text(json.dumps({"candidate_inputs": [file_provenance(p) for p in candidate.iterdir()]}))
    return candidate, audit


def test_consistent_exclusions_and_terminal_stop_removal(tmp_path):
    candidate, audit = data(tmp_path)
    output = tmp_path / "output"
    result = prepare(candidate, output, audit)
    assert result["proteins"] == 1
    assert result["input_exclusions"] == {"OFF": ["off"], "exposed_genus": ["exposed"], "internal_stop": ["b"]}
    assert (output / "input/Vpolyspora.fasta").read_text() == ">a\nAAAA\n"
    assert json.loads((output / "reference_groups.json").read_text()) == {"Pillar00001": ["a"]}


def test_changed_snapshot_rejected(tmp_path):
    candidate, audit = data(tmp_path)
    (candidate / "README").write_text("changed")
    with pytest.raises(ValueError, match="snapshot changed"):
        prepare(candidate, tmp_path / "output", audit)


def test_existing_output_not_overwritten(tmp_path):
    candidate, audit = data(tmp_path)
    output = tmp_path / "output"
    output.mkdir()
    with pytest.raises(ValueError, match="Refusing to overwrite"):
        prepare(candidate, output, audit)
