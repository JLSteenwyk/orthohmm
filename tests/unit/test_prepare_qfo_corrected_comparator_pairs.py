import pytest

from benchmark_tools.prepare_qfo_corrected_comparator_pairs import validate_admission, convert


def admission():
    return {"status": "corrected_proteinortho_native_graph_admitted_for_conversion",
            "accuracy_evaluated": False, "graph_inventory": {"input_species": 78,
                "input_accessions": 984137, "species_sections": 3003, "validated_pair_rows": 1}}


def test_corrected_admission():
    assert validate_admission("proteinortho", admission()) == 1


@pytest.mark.parametrize("key,value", [("input_species", 77), ("input_accessions", 976504),
                                      ("species_sections", 3002), ("validated_pair_rows", 0),
                                      ("validated_pair_rows", True)])
def test_reject_wrong_inventory(key, value):
    report = admission()
    report["graph_inventory"][key] = value
    with pytest.raises(ValueError):
        validate_admission("proteinortho", report)


def test_reject_unadmitted():
    report = admission()
    report["status"] = "running"
    with pytest.raises(ValueError):
        validate_admission("proteinortho", report)


@pytest.mark.parametrize("count", [1, 2])
def test_native_graph_converter_and_count_binding(tmp_path, count):
    graph = tmp_path / "qfo.proteinortho-graph"
    graph.write_text("# s1.fasta\ts2.fasta\nb\ta\t0\t50\t0\t50\n")
    report = admission()
    report["graph_inventory"]["validated_pair_rows"] = count
    report["native_outputs"] = {"pairs": {"path": str(graph)}}
    output = tmp_path / "partial.tsv"
    if count == 1:
        assert convert("proteinortho", report, output) == (1, 0)
        assert output.read_text() == "a\tb\n"
    else:
        with pytest.raises(ValueError, match="count differs"):
            convert("proteinortho", report, output)


def test_sonic_accounting():
    report = {"status": "corrected_sonic_native_tables_admitted_for_conversion", "accuracy_evaluated": False,
              "table_inventory": {"input_species": 78, "input_accessions": 984137,
                "species_pair_tables": 3003, "distinct_pairs": 2, "raw_relations": 3, "duplicate_relations": 1}}
    assert validate_admission("sonic", report) == 2
    report["table_inventory"]["raw_relations"] = 4
    with pytest.raises(ValueError, match="accounting"):
        validate_admission("sonic", report)


def test_sonic_count_binding(monkeypatch, tmp_path):
    from benchmark_tools import prepare_qfo_corrected_comparator_pairs as module
    report = {"status": "corrected_sonic_native_tables_admitted_for_conversion", "accuracy_evaluated": False,
              "pair_directory": str(tmp_path), "table_inventory": {"input_species": 78, "input_accessions": 984137,
                "species_pair_tables": 3003, "distinct_pairs": 2, "raw_relations": 3, "duplicate_relations": 1}}
    monkeypatch.setattr(module, "sonic_pairs", lambda *args: (3003, 2, 0))
    with pytest.raises(ValueError, match="duplicate"):
        convert("sonic", report, tmp_path / "partial.tsv")


@pytest.mark.parametrize("mapped", [True, False])
def test_stage_success_or_mapping_failure_is_preserved(monkeypatch, tmp_path, mapped):
    import gzip
    import json
    from benchmark_tools import prepare_qfo_corrected_comparator_pairs as module
    graph = tmp_path / "qfo.proteinortho-graph"
    graph.write_text("# s1.fasta\ts2.fasta\na\tb\t0\t50\t0\t50\n")
    mapping = tmp_path / "mapping.json.gz"
    with gzip.open(mapping, "wt") as stream:
        json.dump({"mapping": {"a": 1, **({"b": 2} if mapped else {})}}, stream)
    manifest = tmp_path / "admission.json"
    manifest.write_text("{}")
    environment = tmp_path / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment.parent.mkdir(parents=True)
    environment.write_text("{}")
    report = admission()
    report.update(source=module.record(graph), checked_records=[module.record(graph)],
                  native_outputs={"pairs": module.record(graph)})
    monkeypatch.setattr(module, "read_frozen", lambda path, sha: report if path == manifest else
                        {"reference_files": [module.record(mapping)]})
    directory = tmp_path / "benchmarks/results/qfo_corrected_comparator_pairs_v1/proteinortho"
    if mapped:
        result = module.prepare(tmp_path, "proteinortho", manifest, "fixture")
        assert result["status"] == "corrected_comparator_pairs_prepared_unscored"
        assert result["total_pairs"] == result["retained_pairs"] == 1
        assert result["removed_mapping_pairs"] == 0
        assert (directory / "pairs.tsv").read_bytes() == (directory / "pairs.qfo.tsv").read_bytes()
        with pytest.raises(FileExistsError):
            module.prepare(tmp_path, "proteinortho", manifest, "fixture")
    else:
        with pytest.raises(ValueError, match="mapping loss"):
            module.prepare(tmp_path, "proteinortho", manifest, "fixture")
        assert json.loads((directory / "results.json").read_text())["status"] == "failed"
        assert (directory / "pairs.partial.tsv").exists()
        assert not (directory / "pairs.tsv").exists()
