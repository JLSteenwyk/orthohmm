import json

import pytest

from benchmark_tools.derive_simulation_conditions import derive
from benchmark_tools.benchmark_production import file_record
from tests.unit.test_zombi_truth import fixture_run


def test_export_provenance_and_inapplicable_conditions(tmp_path):
    run, output = tmp_path / "run", tmp_path / "derived"
    fixture_run(run)
    report = derive(run, output, 10)
    assert report["accuracy_computed"] is False
    assert report["conditions"]["uneven_taxa"]["status"] == "inapplicable"
    assert report["conditions"]["taxon_count_control"]["status"] == "inapplicable"
    target = output / "missing20"
    truth = json.loads((target / "truth.json").read_text())
    assert truth["ortholog_pair_count"] == 1
    assert truth["extant_genes"] == 2
    assert truth["selection"]["removed_genes"] == []
    assert (target / "input/A.fasta").read_text() == ">F1__A_2\nACDE\n"
    for record in truth["inputs"]:
        assert record == file_record(target / record["path"], target)
    before = (output / "manifest.json").read_bytes()
    with pytest.raises(FileExistsError):
        derive(run, output, 10)
    assert (output / "manifest.json").read_bytes() == before


def test_failed_input_validation_does_not_create_export(tmp_path):
    run, output = tmp_path / "run", tmp_path / "derived"
    fixture_run(run)
    (run / "S/1_complete.fasta").unlink()
    with pytest.raises(ValueError):
        derive(run, output, 10)
    assert not output.exists()
