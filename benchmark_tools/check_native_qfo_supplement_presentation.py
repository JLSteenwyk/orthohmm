"""Check generated tables and decoded PDF figures, without scientific replay."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import fitz

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

TABLES = ("scores", "swiss_intervals", "vgnc_counts", "vgnc_transitions", "swiss_transitions",
          "graph_support", "strata_differences", "functional_overlap")
FLOAT_COLUMNS = dict(scores=(2, 3), swiss_intervals=(1, 2, 3), vgnc_counts=(4, 5, 6), strata_differences=(2, 3, 4))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def text(node):
    if isinstance(node, list):
        return "".join(text(child) for child in node)
    if not isinstance(node, dict):
        return ""
    kind = node.get("t")
    if kind == "Str":
        return node["c"]
    if kind in ("Space", "SoftBreak", "LineBreak"):
        return " "
    if kind == "Code":
        return node["c"][1]
    return text(node.get("c", []))


def table_rows(block):
    require(block["t"] == "Table" and len(block["c"]) == 6, "Unexpected Pandoc table schema")
    rows = []
    for body in block["c"][4]:
        require(not body[2], "Unexpected intermediate header")
        for row in body[3]:
            cells = row[1]
            require(all(cell[2:4] == [1, 1] for cell in cells), "Unexpected spanned cell")
            rows.append([text(cell[4]) for cell in cells])
    return rows


def expected_rows(path, name):
    with path.open(newline="") as stream:
        rows = list(csv.reader(stream, delimiter="\t"))[1:]
    if name == "swiss_transitions":
        rows = [row for row in rows if int(row[2])]
    if name == "strata_differences":
        rows = [row for row in rows if row[0] in ("higher_entropy", "lower_entropy", "short_relative", "not_short_relative")]
    for row in rows:
        for index in FLOAT_COLUMNS.get(name, ()):
            if row[index] != "Unavailable":
                row[index] = f"{float(row[index]):.6f}"
    return rows


def decoded(pixmap):
    if pixmap.alpha:
        pixmap = fitz.Pixmap(pixmap, 0)
    return dict(width=pixmap.width, height=pixmap.height, channels=pixmap.n,
        sha256=hashlib.sha256(pixmap.samples).hexdigest())


def review(assembly_path, assets_path, print_path, pdf_path):
    direct = [record(path) for path in (assembly_path, assets_path, print_path, pdf_path, __file__)]
    assembly, assets, printed = [json.loads(Path(ref["path"]).read_text()) for ref in direct[:3]]
    require(assembly["schema"] == "native_qfo_supplement_assembly_v1" and assembly["publication_ready"] is False
        and assembly["scientific_evidence_replayed"] is False, "Changed assembly scope")
    require(assets["status"] == "manuscript_review_rendered" and printed["status"] == "verified_html_printed"
        and (printed["pdf"]["bytes"], printed["pdf"]["sha256"]) == (direct[3]["bytes"], direct[3]["sha256"]),
        "PDF differs from print receipt")
    checked = [*direct, *assembly["outputs"], *assembly["checked_records"], assets["html"],
        *assets["sources"], *assets["targets"]]
    for ref in checked:
        check(ref)
    manuscript = next(Path(ref["path"]) for ref in assembly["outputs"] if ref["path"].endswith("/supplement.md"))
    ast = json.loads(subprocess.check_output(assets["parse_command"], text=True, timeout=60))
    require(Path(assets["parse_command"][-1]).resolve() == manuscript.resolve(), "Wrong parsed manuscript")
    blocks = [block for block in ast["blocks"] if block["t"] == "Table"]
    require(len(blocks) == len(TABLES), "Wrong displayed table count")
    expected = {name: expected_rows(manuscript.parent / (name + ".tsv"), name) for name in TABLES}
    for name, block in zip(TABLES, blocks):
        require(table_rows(block) == expected[name], "Markdown/TSV mismatch: " + name)
    figure_images = [ref for ref in assets["targets"] if ref["path"].endswith(".png")]
    require(len(figure_images) == 2, "Wrong source figure count")
    with fitz.open(pdf_path) as document:
        page_text = [" ".join(page.get_text().split()) for page in document]
        rows_checked, image_pixels, image_refs = 0, {}, set()
        for name in TABLES:
            for row in expected[name]:
                phrase = " ".join(" ".join(row).split())
                require(any(phrase in page for page in page_text), "Printed table row absent: " + name + " " + phrase)
                rows_checked += 1
        for page in document:
            for item in page.get_images(full=True):
                image_refs.add(item[0])
        for xref in image_refs:
            image_pixels[xref] = decoded(fitz.Pixmap(document, xref))
        matches = []
        for ref in figure_images:
            original = decoded(fitz.Pixmap(ref["path"]))
            matching = [xref for xref, pixels in image_pixels.items() if pixels == original]
            require(len(matching) == 1, "Source figure pixels not exactly present in PDF")
            matches.append(dict(source=ref, pdf_xref=matching[0], decoded_pixels=original))
        pages = len(document)
    for ref in checked:
        check(ref)
    return dict(schema="native_qfo_supplement_presentation_readback_v1", source=record(__file__),
        checked_records=direct, displayed_tables=len(TABLES), displayed_rows=rows_checked,
        table_row_counts={name: len(rows) for name, rows in expected.items()}, pdf_pages=pages,
        exact_decoded_figure_matches=matches, generator_imported=False, scientific_evidence_replayed=False,
        publication_ready=False, visual_review_complete=False,
        limitations=["Checks presentation correspondence to source-bound generated TSVs, not a new raw scientific admission.",
            "Table text and exact embedded image pixels do not certify readability, nonoverlap or biological correctness."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("assembly", "assets", "print", "pdf", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = review(args.assembly, args.assets, vars(args)["print"], args.pdf)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(dict(displayed_tables=result["displayed_tables"], displayed_rows=result["displayed_rows"],
        figure_matches=len(result["exact_decoded_figure_matches"]))))
