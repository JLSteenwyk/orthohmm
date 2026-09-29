import json

import pytest

from benchmark_tools.inspect_gerstein_treefam import inspect, rows
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture(tmp_path):
    paths = []
    for name, gene in (("human", "H1"), ("fly", "F1"), ("worm", "F1")):
        path = tmp_path / f"tf.{name}.genes.raw"
        path.write_text(f"# {name}\n1\t{gene}\ttranscript\t1\tgene\t1\tsymbol\t\n")
        paths.append(path)
    for suffix in ("raw", "all"):
        path = tmp_path / f"tf.human_fly_worm.fam_genes.{suffix}"
        path.write_text("# header\nTF1\tH1\tF1\tmissing\n")
        paths.append(path)
    (tmp_path / "download.json").write_text(json.dumps(dict(artifacts=[dict(file=record(p)) for p in paths])))
    return tmp_path


def test_duplicate_species_contents_reported_not_admitted(tmp_path):
    result = inspect(fixture(tmp_path))
    assert result["fly_and_worm_data_rows_identical"] is True
    assert result["memberships"]["raw"]["ids_absent_from_downloaded_gene_tables"] == 1
    assert result["scientific_use_admitted"] is False


def test_changed_download_rejected(tmp_path):
    directory = fixture(tmp_path)
    (directory / "tf.fly.genes.raw").write_text("changed\n")
    with pytest.raises(ValueError):
        inspect(directory)


def test_tab_parser_preserves_empty_trailing_field(tmp_path):
    path = tmp_path / "data"
    path.write_text("# comment\na\tb\t\n")
    assert rows(path) == [["a", "b", ""]]
