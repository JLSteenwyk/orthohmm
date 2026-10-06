import csv
import json
from pathlib import Path
import subprocess

import pytest

from benchmark_tools import check_native_qfo_supplement_presentation as checker


def test_structured_pandoc_table_readback(tmp_path):
    path = tmp_path / "table.md"
    path.write_text("| Metric | Value |\n|---|---|\n| F1 | -0.046418 |\n| `not_scored` | 16004 |\n")
    ast = json.loads(subprocess.check_output(["pandoc", "--from=markdown", "--to=json", str(path)], text=True))
    table = next(block for block in ast["blocks"] if block["t"] == "Table")
    assert checker.table_rows(table) == [["F1", "-0.046418"], ["not_scored", "16004"]]
    table["c"][4][0][3][0][1][0][2] = 2
    with pytest.raises(ValueError, match="spanned"):
        checker.table_rows(table)


@pytest.mark.parametrize("name,rows,expected", (
    ("scores", [["VGNC", "F1", ".666833533", ".898184579"]], [["VGNC", "F1", "0.666834", "0.898185"]]),
    ("swiss_transitions", [["TP", "FN", 463], ["FN", "TP", 0]], [["TP", "FN", "463"]]),
    ("strata_differences", [["higher_entropy", 9, .1, .3, -.01], ["empty", 0, "Unavailable", "Unavailable", "Unavailable"]],
        [["higher_entropy", "9", "0.100000", "0.300000", "-0.010000"]]),
))
def test_display_filter_and_numeric_serialization(tmp_path, name, rows, expected):
    path = tmp_path / "table.tsv"
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["header"])
        writer.writerows(rows)
    assert checker.expected_rows(path, name) == expected


def test_rgb_decoding_preserves_pixels_with_opaque_alpha(tmp_path):
    import fitz

    document = fitz.open()
    page = document.new_page(width=20, height=20)
    page.draw_rect(page.rect, color=(.2, .5, .8), fill=(.2, .5, .8))
    rgb, rgba = page.get_pixmap(alpha=False), page.get_pixmap(alpha=True)
    assert checker.decoded(rgb) == checker.decoded(rgba)
    assert checker.decoded(rgb)["channels"] == 3


def test_checker_does_not_import_generator():
    import ast

    source = ast.parse(Path(checker.__file__).read_text())
    assert not any(isinstance(node, ast.ImportFrom) and node.module
        and "assemble_native_qfo_supplement" in node.module for node in ast.walk(source))
